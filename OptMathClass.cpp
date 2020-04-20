/******************************************************************************
File     : OptMathClass.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The OptMathClass is used to compute mathematical measures, namely 1st and 2nd 
order derivatives, that are used in certain optimization algorithms.

Version History
08-22-03    lsm   created file
03-05-04    lsm   Using ReadParams()
03-24-04    lsm   Added CalcOptimalStepSize() option
07-08-04    lsm   Switch to using ParameterABC
08-12-04    lsm   Removed ModelBackup class 
                  Added m_pXxxxPoint to reduce memory fragmentation
08-18-04    lsm   Added metrics collection and reporting.
03-21-05    lsm   Added support for parameter-specific relative increments
01-01-07    lsm   Added support for additional FD increment types. OptMathClass
                  now uses abstract model base class (ModelABC).
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "OptMathClass.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObservationGroup.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
CTOR()

Setup pointer to model and use the model's information about parameters to 
size and initialize the various matrices and vectors.
******************************************************************************/
OptMathClass::OptMathClass(ModelABC * pModel)
{
   int i;
   ParameterGroup * pParamGroup;

   m_DiffType  = FD_FORWARD;
   m_MinInc = NEARLY_ZERO;
            
   m_pModel = pModel;
   pParamGroup = m_pModel->GetParamGroupPtr();

   m_NumParams = pParamGroup->GetNumParams();

   NEW_PRINT("double", m_NumParams);
   m_pDiffInc  = new double[m_NumParams];
   MEM_CHECK(m_pDiffInc);

   m_DiffIncType = FD_RANGE_REL;
   for(i = 0; i < m_NumParams; i++){ m_pDiffInc[i] = 0.001;}

   NEW_PRINT("double", m_NumParams);
   m_pGrad  = new double[m_NumParams];

   NEW_PRINT("double", m_NumParams);
   m_pHessPoint = new double[m_NumParams];

   NEW_PRINT("double", m_NumParams);
   m_pGradPoint = new double[m_NumParams];

   NEW_PRINT("double", m_NumParams);
   m_pStepPoint = new double[m_NumParams];

   NEW_PRINT("double", m_NumParams);
   m_pDiffPoint = new double[m_NumParams];

   NEW_PRINT("double *", m_NumParams);
   m_pHess  = new double*[m_NumParams];
   MEM_CHECK(m_pHess);

   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pHess[i] = new double[m_NumParams];
   }
   MEM_CHECK(m_pHess[i-1]);

   //configuration file can override certain defaults
   InitFromFile(GetInFileName());

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
WriteMetrics()

Reports on the setup of the Math Class and also various run-time metrics.
******************************************************************************/
void OptMathClass::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nFinite Difference Metrics\n");
   fprintf(pFile, "Difference Type    : ");

   if     (m_DiffType == FD_FORWARD){ fprintf(pFile, "Forward\n");}   
   else if(m_DiffType == FD_OUT_CEN){ fprintf(pFile, "Outside Central\n");}
   else if(m_DiffType == FD_PAR_CEN){ fprintf(pFile, "Parabolic Central\n");}
   else if(m_DiffType == FD_FIT_CEN){ fprintf(pFile, "Best-fit Central\n");}
   else                             { fprintf(pFile, "Unknown\n");}

   fprintf(pFile, "Increment Type    : ");

   if     (m_DiffIncType == FD_RANGE_REL){ fprintf(pFile, "Range-Relative\n");}   
   else if(m_DiffIncType == FD_VALUE_REL){ fprintf(pFile, "Value-Relative\n");}
   else if(m_DiffIncType == FD_ABSOLUTE) { fprintf(pFile, "Absolute\n");}
   else if(m_DiffIncType == FD_OPTIMAL)  { fprintf(pFile, "Optimal\n");}
   else                                  { fprintf(pFile, "Range-Relative\n");}
 
   fprintf(pFile, "Finite Difference Increments\n");
   for(int i = 0; i < m_NumParams; i++)
   {
      fprintf(pFile, "%-12s : ", GetParameterName(i));
      if(m_DiffIncType != FD_OPTIMAL){
         fprintf(pFile, "%lf\n", m_pDiffInc[i]);}
      else{
         fprintf(pFile, "optimal\n");}
   }
   fprintf(pFile, "Finite Difference Mimumum Increment : %E\n", m_MinInc);

   fprintf(pFile, "Hessian Evals      : %d\n", m_HessCount);
   fprintf(pFile, "Gradient Evals     : %d\n", m_GradCount);
   fprintf(pFile, "Derivative Evals   : %d\n", m_DiffCount);
   fprintf(pFile, "Optimal Step Evals : %d\n", m_StepCount);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void OptMathClass::InitFromFile(IroncladString pMathFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char * pTmp;
   int i,j;

   m_DiffCount = 0;
   m_StepCount = 0;
   m_GradCount = 0;
   m_HessCount = 0;

   //read in math configuration
   pFile = fopen(pMathFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open math config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginMathAndStats", pMathFileName) == true)
   {
      FindToken(pFile, "EndMathAndStats", pMathFileName);
      rewind(pFile);

      FindToken(pFile, "BeginMathAndStats", pMathFileName);
      line = GetNxtDataLine(pFile, pMathFileName);
      while(strstr(line, "EndMathAndStats") == NULL)
      {
         if(strstr(line, "DiffType") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            strcpy(line, tmp2);
            MyStrLwr(line);
            if(strstr(line, "forward") != NULL) {m_DiffType = FD_FORWARD;}
            else if(strstr(line, "outside") != NULL) {m_DiffType = FD_OUT_CEN;}
            else if(strstr(line, "parabolic") != NULL) {m_DiffType = FD_PAR_CEN;}
            else if(strstr(line, "best-fit") != NULL) {m_DiffType = FD_FIT_CEN;}
         }/*end if() */
         else if(strstr(line, "DiffIncType") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            strcpy(line, tmp2);
            MyStrLwr(line);
            if     (strstr(line, "range-relative") != NULL) {m_DiffIncType = FD_RANGE_REL;}
            else if(strstr(line, "value-relative") != NULL) {m_DiffIncType = FD_VALUE_REL;}
            else if(strstr(line, "absolute")       != NULL) {m_DiffIncType = FD_ABSOLUTE;}
            else if(strstr(line, "optimal")        != NULL) {m_DiffIncType = FD_OPTIMAL;}
         }/*end if() */
         else if(strstr(line, "DiffRelIncrement") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("DiffRelIncrement");

            for(i = 0; i < m_NumParams; i++)
            {
                j = ExtractString(pTmp, tmp);                
                m_pDiffInc[i] = atof(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            //if only parial list, apply first value to remaining params
            for(i = i; i < m_NumParams; i++)
            { 
               m_pDiffInc[i] = m_pDiffInc[0];
            }

            /* this keyword (DiffRelIncrement) is range-relative */
            m_DiffIncType = FD_RANGE_REL;
         }/*end else if() */         
         else if(strstr(line, "DiffIncrement") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("DiffIncrement");

            for(i = 0; i < m_NumParams; i++)
            {
                j = ExtractString(pTmp, tmp);
                m_pDiffInc[i] = atof(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            //if only parial list, apply first value to remaining params
            for(i = i; i < m_NumParams; i++)
            { 
               m_pDiffInc[i] = m_pDiffInc[0];
            }
         }/*end else if() */
         else if(strstr(line, "DiffMinIncrement") != NULL)
         {
            sscanf(line, "DiffMinIncrement %lf", &m_MinInc);
         }/*end else if() */
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pMathFileName);
      } /* end while() */
   }/* end if() */   

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
Destroy()

Free up arrays.
******************************************************************************/
void OptMathClass::Destroy(void)
{
   int i;   

   delete [] m_pGrad;   

   for(i = 0; i < m_NumParams; i++){delete [] m_pHess[i];}
   delete [] m_pHess;

   delete [] m_pHessPoint;
   delete [] m_pStepPoint;
   delete [] m_pGradPoint;
   delete [] m_pDiffPoint;

   delete [] m_pDiffInc;
   
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
CalcHessian()

Calculate the Hessian matrix (2nd order partial derivatives). This method uses
forward differences for the 2nd order approximation. It calculates the 
Hessian at the point of the most recent model run. After completion of the
Hessian calculation, the model is rerun at the initial location to ensure that
the system remains in a consistent state.
******************************************************************************/
Unchangeable2DArray OptMathClass::CalcHessian(void)
{
   int *iMap, * jMap;
   double * pFij, * pFi, * Fall, * pdx;
   double F, Fij, Fi, Fj, next;
   double *x = m_pHessPoint, dxi, dxj;
   int i, j, ij, k, num_evals, id, np;
   ParameterGroup * pParamGroup;

   x = new double[m_NumParams];
   pdx = new double[m_NumParams];

   //initialize point
   pParamGroup = m_pModel->GetParamGroupPtr();
   pParamGroup->ReadParams(x);

   //assign delta x values, reversing direction if near boundaries
   for(i = 0; i < m_NumParams; i++)
   {      
      if(m_DiffIncType == FD_ABSOLUTE)
      {
         dxi = fabs(m_pDiffInc[i]);
      }
      else if(m_DiffIncType == FD_RANGE_REL)
      {
         double range = pParamGroup->GetParamPtr(i)->GetUprBnd() - pParamGroup->GetParamPtr(i)->GetLwrBnd();
         dxi = fabs(m_pDiffInc[i])*range;
      }
      else //value-relative
      {
         dxi = fabs(m_pDiffInc[i])*x[i];
         if(dxi < m_MinInc) dxi = m_MinInc;
      }

      //trick from NR in C
      next = x[i] + dxi;
      dxi = next - x[i];
      pdx[i] = dxi;
      if((x[i]+2*dxi) > pParamGroup->GetParamPtr(i)->GetUprBnd())
            pdx[i] = -dxi;      
   }/* end for() */

   //create storage for list of required objective function values
   num_evals = 1+m_NumParams+((m_NumParams*(m_NumParams+1))/2);
   Fall = new double[num_evals];
   iMap = new int[num_evals];
   jMap = new int[num_evals];
   for(i = 0; i < num_evals; i++)
   {
      Fall[i] = 0.00;
      iMap[i] = 0;
      jMap[i] = 0;
   }/* end for() */
   pFi = &(Fall[0]); //single perturbations
   pFij = &(Fall[m_NumParams]); //double perterbations

   //assign index mappings
   k = m_NumParams;
   for(i = 0; i < m_NumParams; i++)
   {
      for(j = i; j < m_NumParams; j++)
      {
         iMap[k] = i;
         jMap[k] = j;
         k++;
      }/* end for() */
   }/* end for() */

   //compute objective function values, possibly in parallel
   MPI_Comm_size(MPI_COMM_WORLD, &np);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   for(i = id; i < num_evals; i+=np)
   {
      if(i < m_NumParams)
      {
         dxi = pdx[i];
         x[i] += dxi;
         pParamGroup->WriteParams(x);
         pFi[i] = m_pModel->Execute();
         x[i] -= dxi;
      }/* end if() */
      else if (i == (num_evals-1))
      { 
        pParamGroup->WriteParams(x);
        F = m_pModel->Execute();
        Fall[i] = F;
      }
      else
      {   
         dxi = pdx[iMap[i]];
         dxj = pdx[jMap[i]];
         x[iMap[i]] += dxi;
         x[jMap[i]] += dxj;
         pParamGroup->WriteParams(x);
         ij = i - m_NumParams;
         pFij[ij] = m_pModel->Execute();
         x[iMap[i]] -= dxi;
         x[jMap[i]] -= dxj;
      }/* end else() */
   }/* end for() */

   //if parallel, collect results
   if(np > 1)
   {
      double * pTmp = new double[num_evals];
      for(i = 0; i < num_evals; i++){ pTmp[i] = 0.00;}
      MPI_Reduce(Fall, pTmp, num_evals,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      for(i = 0; i < num_evals; i++){ Fall[i] = pTmp[i];}
      delete [] pTmp;
      F = Fall[num_evals - 1];
   }/* end if() */

   //compute Hessian matrix
   k = 0;
   for(i = 0; i < m_NumParams; i++)
   {
      for(j = i; j < m_NumParams; j++)
      {
         dxi = pdx[i];
         dxj = pdx[j];
         Fij = pFij[k];
         k++;
         Fj = pFi[j];
         Fi = pFi[i];
         m_pHess[i][j] = (Fij - Fi - Fj + F) / (dxi * dxj);
         //Hessian is symmetric
         if(i != j) m_pHess[j][i] =  m_pHess[i][j];         
      }/* end for() */
   }/* end for() */

   delete [] x;
   delete [] pdx;
   delete [] Fall;
   delete [] iMap;
   delete [] jMap;
      
   return m_pHess;
}/* end CalcHessian() */

/******************************************************************************
CalcDerivative()

Calculates the derivative with respect to the given paramter (indicated by the 
parmIdx argument) using finite differences. The derivative is computed relative
to the point at which the model was last executed. 

NOTE: After completion of the derivative calculation, the model is NOT rerun 
at the initial location, therefore this routine will leave the system in an 
inconsistent state. Routines that call CalcDerivative() must handle the 
restoration of the system.

If a better minimum than fmin is found, fmin and pmin are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
double OptMathClass::CalcDerivative(int parmIdx, double * fmin, double * pmin)
{
   int j;
   double midParm, lhsParm, rhsParm;
   double midObj, lhsObj, rhsObj;
   double dx, diff, upr, lwr, cur, next;
   double Sxy, Sx, Sy, Sxx; //central diff. best-fit parameters   
   FiniteDiffType dType;
   FiniteDiffIncType dIncType;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   dType = m_DiffType;
   dIncType = m_DiffIncType;

//if FD calculation is ~0.00, retry using an alternative increment type
retry:
   //initialize current state of model      
   pGroup = m_pModel->GetParamGroupPtr();
   //read in current point
   pGroup->ReadParams(m_pDiffPoint);
   //assign middle obj. function 
   midObj = m_pModel->GetObjFuncVal();

   /*----------------------------------- 
   compute left-hand side and right-hand 
   side locations at which the objective
   function will be calculated.
   ------------------------------------*/   
   //retrieve parameter limits
   pParam = pGroup->GetParamPtr(parmIdx);
   cur = pParam->GetEstVal();
   upr = pParam->GetUprBnd();
   lwr = pParam->GetLwrBnd();
   //setup step size
   midParm = m_pDiffPoint[parmIdx];

   j = parmIdx;
   if(dIncType == FD_OPTIMAL)
   {
      dx = CalcOptimalStepSize(j);
   }
   else if(dIncType == FD_RANGE_REL)
   {
      dx = fabs(m_pDiffInc[j]*(upr-lwr));
   }
   else if(dIncType == FD_VALUE_REL)
   {
      dx = MyMax(fabs(m_pDiffInc[j]*cur), m_MinInc);
   }
   else if(dIncType == FD_ABSOLUTE)
   {
      dx = fabs(m_pDiffInc[j]);
   }
   else //default to range-relative increment
   {
      dx = fabs(m_pDiffInc[j]*(upr-lwr));
   }
   //trick from NR in C   
   next = cur + dx;
   dx = next - cur;
   
   //set perterbation steps
   if(dType == FD_FORWARD)
   {
      rhsParm =  midParm + dx;
      lhsParm = midParm - dx;
   }
   else //(dType != FD_FORWARD)
   //only take half steps if not using forward diff.
   {
      rhsParm =  midParm + (0.5 * dx);
      lhsParm = midParm - (0.5 * dx);
   }/* end if() */

   //avoid exceeding parameter limits
   if(rhsParm > upr) 
   { 
      //switch direction and difference type
      dType = FD_FORWARD;
      dx *= -1;
      rhsParm = midParm + dx;
   } /* end if() */
   if(lhsParm < lwr)
   { 
      //switch difference type
      dType = FD_FORWARD;
      rhsParm = midParm + dx;
   }/* end if() */
   
   //compute rhs obj. func.
   m_pDiffPoint[parmIdx] = rhsParm;
   pGroup->WriteParams(m_pDiffPoint);
   rhsObj = m_pModel->Execute();
   /* update optimal, if appropriate */
   if((fmin != NULL) && (rhsObj < *fmin))
   {
      *fmin = rhsObj;
      pGroup->ReadParams(pmin);
   }
   m_DiffCount++;
   
   //compute lhs obj. func., if needed
   if(dType != FD_FORWARD)
   {
      m_pDiffPoint[parmIdx] = lhsParm;
      pGroup->WriteParams(m_pDiffPoint);
      lhsObj = m_pModel->Execute();
      /* update optimal, if appropriate */
      if((fmin != NULL) && (lhsObj < *fmin))
      {
         *fmin = lhsObj;
         pGroup->ReadParams(pmin);
      }
      m_DiffCount++;
   }/* end if() */

   /*---------------------------------------
   Compute the partial derivative.
   ----------------------------------------*/
   switch(dType)
   {
      case(FD_OUT_CEN) : //outside central
         dx = (rhsParm - lhsParm);
         diff = ((rhsObj - lhsObj) / dx);
         break;
      case(FD_PAR_CEN) : //parabolic central
         dx = (rhsParm - lhsParm);
         double x1,x2,x3,F1,F2,F3,denom;
         double a, b, c, dF1, dF2, dx1, dx2;

         x1 = lhsParm; x2 = midParm; x3 = rhsParm;
         F1 = lhsObj;  F2 = midObj;  F3 = rhsObj;

         denom = (((x3*x3 - x1*x1))-((x2+x1)*(x3-x1)));
         dF1 = F3-F1;
         dF2 = F2-F1;
         dx1 = x3-x1;
         dx2 = x2-x1;

         a = (dF1-(dF2*dx1)/dx2)/denom;
         b = (dF2-a*((x2*x2)-(x1*x1)))/dx2;
         c = F1-a*x1*x1-b*x1;

         diff = 2*a*midParm + b;

         if(fabs(diff) < NEARLY_ZERO)
         {
            if(dIncType != FD_ABSOLUTE) //try switching increment type
            {
               dIncType = FD_ABSOLUTE;

               //semi-restore the model (for next time around)   
               m_pDiffPoint[parmIdx] = midParm;
               pGroup->WriteParams(m_pDiffPoint);   
               m_pModel->SetObjFuncVal(midObj);

               goto retry;
            }
            else //use outside central
            {
               diff = ((rhsObj - lhsObj) / dx);
            }
         }
         break;
      case(FD_FIT_CEN) : //Least-Squares best-fit central
         /*-------------------------------------------------------               
         The derivative of (y = bx + a) using Least-Sqaures is: 
            dy/dx = b = (SSxy - SxSy)/(SSxx - (Sx)^2)
            where, 
            S = #points (3),
            Sxy = (x1y1 + x2y2 + x3y3),
            Sx = (x1 + x2 + x3),
            Sy = (y1 + y2 + y3),
            Sxx = (x1^2 + x2^2 + x3^2),
            x1,x2,x3 -> parameters, 
            y1,y2,y3 -> objective function values
            and equal variance has been given to each data point.
         -------------------------------------------------------*/
         Sxy = ((lhsObj*lhsParm)+(midObj*midParm)+(rhsObj*rhsParm));
         Sx = (lhsParm + midParm + rhsParm);
         Sy = (lhsObj + midObj + rhsObj);
         Sxx=(lhsParm*lhsParm)+(midParm*midParm)+(rhsParm*rhsParm);
         diff = ((3.00 * Sxy) - (Sx * Sy)) / ((3.00 * Sxx) - (Sx * Sx));
         break;
      case(FD_FORWARD) : //forward difference
      default:
         diff = ((rhsObj - midObj) / dx);
         break;         
   }/* end switch() */     

   //semi-restore the model (for next time around)   
   m_pDiffPoint[parmIdx] = midParm;
   pGroup->WriteParams(m_pDiffPoint);   
   m_pModel->SetObjFuncVal(midObj);

   //if FD calculation is ~0.00, retry using an alternative increment type
   if((fabs(diff) <= NEARLY_ZERO) && (dIncType != FD_RANGE_REL) && (GetProgramType() != GRADIENT_PROGRAM))
   {
      dType = FD_FORWARD;
      dIncType = FD_RANGE_REL;

      goto retry;
   }
  
   return diff;
}/* end CalcDerivative() */

/******************************************************************************
CalcOptimalStepSize()

Calculates the optimal step size using the equations (4) and (5) from Yager, 2004.
"Effects of Model Sensitivty and Nonlinearity on Nonlinear Regression of Ground-
Water Flow".

NOTE: After completion of the optimal step size calculation, the model is NOT 
rerun at the initial location, therefore this routine will leave the system in 
an inconsistent state. Routines that call CalcOptimalStepSize() must handle the 
restoration of the system.
******************************************************************************/
double OptMathClass::CalcOptimalStepSize(int idx)
{
   double delta = 1.00;
   double db, old_db, eps, tmp;
   double Fmid, Fupr, Flwr;
   double bMid;
   double Sjj;   
   
   //read in the point at which the optimal step size will be calculated
   m_pModel->GetParamGroupPtr()->ReadParams(m_pStepPoint);
   bMid = m_pStepPoint[idx];
   Fmid = m_pModel->GetObjFuncVal();

   /*------------------------------------------------------------
   This loop iterates on the step size db, until Yager's criteria
   for optimal step size is met to 3-decimal accuracy.
   -------------------------------------------------------------*/
   eps = 0.001;
   db = old_db = 2.00*sqrt(eps)*fabs(bMid);
   while(delta > eps)
   {
      //forward step
      m_pStepPoint[idx] = bMid + db;
      m_pModel->GetParamGroupPtr()->WriteParams(m_pStepPoint);
      Fupr = m_pModel->Execute();
      m_StepCount++;

      //backward step
      m_pStepPoint[idx] = bMid - db;
      m_pModel->GetParamGroupPtr()->WriteParams(m_pStepPoint);
      Flwr = m_pModel->Execute();
      m_StepCount++;

      /*-----------------------------------------------------
      Revise step size according to Yager's formula, extra
      checks are performed to avoid numerical instability 
      (i.e. taking sqrt() of negative number).
      -------------------------------------------------------*/
      Sjj = (Fupr - 2.00*Fmid + Flwr) / (db*db);
      if(Sjj == 0.00){ db = 2.00*sqrt(eps)*fabs(bMid); break;}
      tmp = (4.00*eps*Fmid)/Sjj;
      if(tmp <= 0.00){ db = 2.00*sqrt(eps)*fabs(bMid); break;}
      db = sqrt(fabs(tmp));      
      delta = fabs(db-old_db);
      old_db = db;
   }/* end if() */

   //semi-restore the model (for next time around)
   m_pStepPoint[idx] = bMid;
   m_pModel->GetParamGroupPtr()->WriteParams(m_pStepPoint);
   m_pModel->SetObjFuncVal(Fmid);

   return db;
}/* end CalcOptimalStepSize() */

/******************************************************************************
CalcGradient()

Cacluates the gradient, relative to the most recent model parameters. After 
completion of the gradient calculation, the model is rerun at the initial 
location to ensure that the system remains in a consistent state.

If a better minimum than fmin is found, fmin and pmin are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
Unchangeable1DArray OptMathClass::CalcGradient(double * fmin, double * pmin)
{
   int i, id, np;
   ParameterGroup * pParamGroup;
   double Finit, Fcur;
 
   //save the design point at which the gradient is to be calculated
   pParamGroup = m_pModel->GetParamGroupPtr();
   pParamGroup->ReadParams(m_pGradPoint);
   Finit = m_pModel->GetObjFuncVal();

   //compute partial derivatives, filling gradient matrix
   MPI_Comm_size(MPI_COMM_WORLD, &np);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   for(i = 0; i < m_NumParams; i++){ m_pGrad[i] = 0.00;}
   for(i = id; i < m_NumParams; i+=np) { m_pGrad[i] = CalcDerivative(i, fmin, pmin);}
   //if parallel, collect results
   if(np > 1)
   {
      double * pTmp = new double[m_NumParams];
      for(i = 0; i < m_NumParams; i++){ pTmp[i] = 0.00;}
      MPI_Reduce(m_pGrad, pTmp, m_NumParams,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      for(i = 0; i < m_NumParams; i++){ m_pGrad[i] = pTmp[i];}
      delete [] pTmp;
   }/* end if() */

   //restore model consistency, unless program type is Gradient_Program
   if(GetProgramType() != GRADIENT_PROGRAM)
   {
      pParamGroup->WriteParams(m_pGradPoint);
      Fcur = m_pModel->Execute();
      m_GradCount++;

      //check model consistency
      if(Fcur != Finit)
      {
         LogError(ERR_MODL_EXE, "CalcGradient() caused model to be inconsistent");
      }   
   }
   
   return m_pGrad;
}/* end CalcGradient() */

/******************************************************************************
Hessian_Program()

Compute the Hessian of the parameter set defined in the input file.
******************************************************************************/
void Hessian_Program(int argc, StringType argv[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("OptMathClass", 1);
   OptMathClass * math = new OptMathClass(model);
   MEM_CHECK(math);

   FILE * pFile;
   int id, j;
   char outName[DEF_STR_SZ];
   char tmp[DEF_STR_SZ];
   const char * inFile = GetOstFileName();
   char * line, * pTok;

   //allocate space for the parameter list
   int num = model->GetParamGroupPtr()->GetNumParams();
   double * pVals = new double[num];

   /* read in user-specified parameter set */
   pFile = fopen(inFile, "r");
   FindToken(pFile, "BeginInitParams", inFile);
   line = GetNxtDataLine(pFile, inFile);

   pTok = line;
   //extract values, one-by-one, making any necessary conversions
   for(int k = 0; k < num; k++)
   {
      j = ExtractString(pTok, tmp);
      j = ValidateExtraction(j, k, num, "Hessian_Program()");
      pTok += j;            
      pVals[k] = model->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
   }/* end for() */
   model->GetParamGroupPtr()->WriteParams(pVals);
   delete [] pVals;

   FindToken(pFile, "EndInitParams", inFile);
   fclose(pFile);

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(outName, "OstOutput%d.txt", id);

   //write setup to file
   pFile = fopen(outName, "w");
   fprintf(pFile, "Ostrich Setup\n");
   fprintf(pFile, "Model: %s\n", model->GetModelStr());
   fprintf(pFile, "Algorithm: Hessian Calculation\n");
   fprintf(pFile, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(pFile, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(pFile, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(pFile, "0\n");}
   else {fprintf(pFile, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fprintf(pFile, "Hessian matrix written to OstHessian.txt\n");
   fclose(pFile);

   //write setup to stdout
   fprintf(stdout, "Ostrich Setup\n");
   fprintf(stdout, "Model: %s\n", model->GetModelStr());
   fprintf(stdout, "Algorithm: Hessian Calculation\n");
   fprintf(stdout, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(stdout, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(stdout, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(stdout, "0\n");}
   else {fprintf(stdout, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fprintf(stdout, "Hessian matrix written to OstHessian.txt\n");

   //model->Execute();
   Unchangeable2DArray pH = math->CalcHessian(); //compute Hessian, possibly in parallel
      
   if(id == 0)
   {
      FILE * pOut = fopen("OstHessian.txt", "w");
      for(int i = 0; i < model->GetParamGroupPtr()->GetNumParams(); i++)
      {
	     for(int j = 0; j < model->GetParamGroupPtr()->GetNumParams(); j++)
		  {
           fprintf(pOut, "%.14E ", pH[i][j]);
        }/* end for() */

        fprintf(pOut, "\n");
      }/* end for() */
      fclose(pOut);
   }/* end if() */
      
   pFile = fopen(outName, "a");
   math->WriteMetrics(pFile);
   fclose(pFile);
   math->WriteMetrics(stdout);

   delete math;
   delete model;
}/* end Hessian_Program() */

/******************************************************************************
Gradient_Program()

Compute the gradient of the parameter set defined in the input file.
******************************************************************************/
void Gradient_Program(int argc, StringType argv[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("OptMathClass", 1);
   OptMathClass * math = new OptMathClass(model);
   MEM_CHECK(math);

   FILE * pFile;
   int id, j;
   char outName[DEF_STR_SZ];
   char tmp[DEF_STR_SZ];
   const char * inFile = GetOstFileName();
   char * line, * pTok;

   //allocate space for the parameter list
   int num = model->GetParamGroupPtr()->GetNumParams();
   double * pVals = new double[num];

   /* read in user-specified parameter set */
   pFile = fopen(inFile, "r");
   FindToken(pFile, "BeginInitParams", inFile);
   line = GetNxtDataLine(pFile, inFile);

   pTok = line;
   //extract values, one-by-one, making any necessary conversions
   for(int k = 0; k < num; k++)
   {
      j = ExtractString(pTok, tmp);
      j = ValidateExtraction(j, k, num, "Gradient_Program()");
      pTok += j;            
      pVals[k] = model->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
   }/* end for() */
   model->GetParamGroupPtr()->WriteParams(pVals);
   delete [] pVals;

   FindToken(pFile, "EndInitParams", inFile);
   fclose(pFile);

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(outName, "OstOutput%d.txt", id);

   //write setup to file
   pFile = fopen(outName, "w");
   fprintf(pFile, "Ostrich Setup\n");
   fprintf(pFile, "Model: %s\n", model->GetModelStr());
   fprintf(pFile, "Algorithm: Gradient Calculation\n");
   fprintf(pFile, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(pFile, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(pFile, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(pFile, "0\n");}
   else {fprintf(pFile, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fprintf(pFile, "Gradient vector written to OstGradient.txt\n");
   fclose(pFile);

   //write setup to stdout
   fprintf(stdout, "Ostrich Setup\n");
   fprintf(stdout, "Model: %s\n", model->GetModelStr());
   fprintf(stdout, "Algorithm: Gradient Calculation\n");
   fprintf(stdout, "Objective Function: %s\n", model->GetObjFuncStr());
   fprintf(stdout, "Number of Parameters: %d\n", model->GetParamGroupPtr()->GetNumParams());
   fprintf(stdout, "Number of Observations: ");
   if(model->GetObsGroupPtr() == NULL){fprintf(stdout, "0\n");}
   else {fprintf(stdout, "%d\n", model->GetObsGroupPtr()->GetNumObs());}
   fprintf(stdout, "Gradient vector written to OstGradient.txt\n");

   model->Execute();
   Unchangeable1DArray pG = math->CalcGradient(NULL, NULL); //compute gradient, possibly in parallel
   
   if(id == 0)
   {
      FILE * pOut = fopen("OstGradient.txt", "w");
      for(int i = 0; i < model->GetParamGroupPtr()->GetNumParams(); i++)
      {
         fprintf(pOut, "%.14E\n", pG[i]);
      }/* end for() */
      fclose(pOut);
   }/* end if() */
   
   pFile = fopen(outName, "a");
   math->WriteMetrics(pFile);
   fclose(pFile);
   math->WriteMetrics(stdout);

   delete math;
   delete model;
}/* end Gradient_Program() */
