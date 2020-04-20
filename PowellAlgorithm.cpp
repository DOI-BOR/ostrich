/******************************************************************************
File      : PowellAlgorithm.cpp
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of Powell's optimization algorithm.

Version History
08-26-03    lsm   created 
02-17-04    lsm   switched implementation
07-08-04    lsm   WriteSetup() is first action of algorithm.
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
08-02-05    lsm   Scrapped NRC implementation, went back to Vanderplaats version.
                  Added status reporting (support for grid computing)
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel. 
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "PowellAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "OptSearchClass.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void PowellAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Initializes parameters, reading user-specified input, if available.
******************************************************************************/
PowellAlgorithm::PowellAlgorithm(ModelABC * pModel)
{
   int i;
   ParameterGroup * pGroup;
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();
   m_pStats = NULL;

   //init. everything to reasonable defaults
   m_MaxIter = 20;
   m_ConvVal = 1E-6; 
   m_AlgCount = 0;
   m_NumRestarts = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;

   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   m_NumDirs = pGroup->GetNumParams();   

   NEW_PRINT("double *", m_NumDirs);
   m_pSearchDirs = new double*[m_NumDirs];
   MEM_CHECK(m_pSearchDirs);

   for(i = 0; i < m_NumDirs; i++)
   { 
      NEW_PRINT("double", m_NumDirs);
      m_pSearchDirs[i] = new double[m_NumDirs];
   }
   MEM_CHECK(m_pSearchDirs[i-1]);

   NEW_PRINT("OptSearchClass", 1);
   m_pSearch = new OptSearchClass(pModel);
   MEM_CHECK(m_pSearch);

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("PowellAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginPowellAlg", pFileName) == true)
   {
      FindToken(inFile, "EndPowellAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginPowellAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndPowellAlg") == NULL)
      {
         if(strstr(line, "ConvergenceVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_ConvVal);
         }
         else if(strstr(line, "MaxIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxIter);
         }

         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default algorithm setup.");
   }/* end else() */

   fclose(inFile);

   IncCtorCount();
}/* end default CTOR */

/******************************************************************************
Destroy()

Frees up memory used by the algorithm.
******************************************************************************/
void PowellAlgorithm::Destroy(void)
{
   int i;

   for(i = 0; i < m_NumDirs; i++){ delete [] m_pSearchDirs[i];}
   delete [] m_pSearchDirs;

   delete m_pStats;
   delete m_pSearch;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Optimize()

Optimize the objective function using the well-known Powell's algorithm.
******************************************************************************/
void PowellAlgorithm::Optimize(void)
{
   StatusStruct pStatus;
   int iter, i, j, n, count;
   ParameterGroup * pGroup;
   double * X; //parameter vector
   double ** S; //search matrix
   double Ftol, Fcur, xmin, Fold, max, tmp, Ftmp;
   double *upr, *lwr, * Xold, * Scur;
  
   count = 0;
   n = m_NumDirs;
   S = m_pSearchDirs;

   double fmin;
   double * pmin;
   NEW_PRINT("double", n);
   pmin = new double[n];

   NEW_PRINT("double", n);
   X = new double[n];

   NEW_PRINT("double", n);
   Xold = new double[n];

   NEW_PRINT("double", n);
   Scur = new double[n];

   NEW_PRINT("double", n);
   upr = new double[n];

   NEW_PRINT("double", n);
   lwr = new double[n];

   Ftol = m_ConvVal;

   MEM_CHECK(lwr);

   //write setup
   WriteSetup(m_pModel, "Powell's Method");

   //initialize search directions
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         if(i == j){S[i][j] = 1.00;}
         else {S[i][j] = 0.00;}
      }/* end for() */
   }/* end for() */

   //initialize X
   pGroup = m_pModel->GetParamGroupPtr();
   pGroup->ReadParams(X);

   //handle warm start
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
      pGroup->ReadParams(X);
   }

   //get uppper and lower bounds
   for(j = 0; j < n; j++)
   {
      upr[j] = pGroup->GetParamPtr(j)->GetUprBnd();
      lwr[j] = pGroup->GetParamPtr(j)->GetLwrBnd();
   }/* end for() */

   //initialize fret   
   Fcur = Fold = m_pModel->Execute();
   m_AlgCount++;

   //write banner and initial result
   WriteBanner(m_pModel, "iter  obj. function  ", "dObjFunc");
   WriteRecord(m_pModel, 0, Fcur, Fcur);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxIter;
   pStatus.pct = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //for(j = 0; j < n; j++){ ptt[j] = pt[j] = p[j];}

   //main loop
   for(iter = 0; iter < m_MaxIter; iter++)
   {
      if(IsQuit() == true){ break;}

      pStatus.curIter = m_CurIter = iter+1;
      
      //loop over all directions
      for(i = 0; i < n; i++)
      {
         //current direction
         for(j = 0; j < n; j++){ Scur[j] = S[i][j];} 

         /*------------------------------------------------------
         normalize search direction so that largest value is 1.00
         -------------------------------------------------------*/
         //first find max value (in absolute terms)
         max = 0.00;
         for(j = 0; j < n; j++)
         { 
            tmp = fabs(Scur[j]);
            if(tmp > max) { max = tmp;}
         }/* end for() */
         //now divide through by max, or reset if max==0
         if(fabs(max) > NEARLY_ZERO)
            for(j = 0; j < n; j++) { Scur[j] /= max;}
         else
         {
            for(j = 0; j < n; j++) 
            { 
               Scur[j] = 0.00;

               for(int k = 0; k < n; k++)
               {
                  if(j == k) S[j][k] = 1.00;
                  else       S[j][k] = 0.00;
               }
            }
            Scur[i] = 1.00;
         }

         //line minimization
		 fmin = Fcur;
         xmin = m_pSearch->CalcStepSize(Scur, &fmin, pmin);
		 if(fmin < Fcur) //found new minimum when computing step size
		 {
            pGroup->WriteParams(pmin);
	        pGroup->ReadParams(X);
		    m_pModel->SetObjFuncVal(fmin);
	        Fcur = fmin;
		 }/* end if() */

         for(j= 0; j < n; j++)
         {
            Scur[j] *= xmin;            
            Xold[j] = X[j];
            X[j] += Scur[j];
            //if it's out of bounds, move half the distance to upr/lwr
            if(X[j] > upr[j]){X[j] = (upr[j]+Xold[j])/2.00; m_NumUprViols++;}            
            if(X[j] < lwr[j]){X[j] = (Xold[j]+lwr[j])/2.00; m_NumLwrViols++;}
         }/* end for() */
         pGroup->WriteParams(X);
         Ftmp = m_pModel->Execute();
         //revert if move didn't improve objective
         if(Ftmp <= Fcur)
         { 
            Fcur = Ftmp;
            //save the revised direction
            for(j = 0; j < n; j++){ S[i][j] = Scur[j];}            
         }
         else
         {
            for(j= 0; j < n; j++){ X[j] = Xold[j];}
            pGroup->WriteParams(X);
            m_pModel->SetObjFuncVal(Fcur);
         }
         m_AlgCount++;
         //end line minimization
      } /* end for() */

      //construct conjugate direction
      for(j = 0; j < n; j++)
      {
         Scur[j] = 0.00;
         for(i = 0; i < n; i++){ Scur[j] += S[i][j];}
      }/* end for() */

      /*------------------------------------------------------
      normalize search direction so that largest value is 1.00
      -------------------------------------------------------*/
      //first find max value (in absolute terms)
      max = 0.00;
      for(j = 0; j < n; j++)
      { 
         tmp = fabs(Scur[j]);
         if(tmp > max) { max = tmp;}
      }/* end for() */
      //now divide through by max, if not zero
      if(fabs(max) > NEARLY_ZERO)
         for(j = 0; j < n; j++) { Scur[j] /= max;}
//      else
//         printf("empty conjugate dir!\n");

      //line minimization
      fmin = Fcur;
      xmin = m_pSearch->CalcStepSize(Scur, &fmin, pmin);
      if(fmin < Fcur) //found new minimum when computing step size
	  {
         pGroup->WriteParams(pmin);
	     pGroup->ReadParams(X);
		 m_pModel->SetObjFuncVal(fmin);
	     Fcur = fmin;
      }/* end if() */

      for(j= 0; j < n; j++)
      {
         Scur[j] *= xmin;            
         Xold[j] = X[j];
         X[j] += Scur[j];
         //if it's out of bounds, move half the distance to upr/lwr
         if(X[j] > upr[j]){X[j] = (upr[j]+Xold[j])/2.00; m_NumUprViols++;}            
         if(X[j] < lwr[j]){X[j] = (Xold[j]+lwr[j])/2.00; m_NumLwrViols++;}
      }/* end for() */
      pGroup->WriteParams(X);
      Ftmp = m_pModel->Execute();
      //revert if move didn't improve objective
      if(Ftmp <= Fcur)
      { 
         Fcur = Ftmp;
         //shift direction matrix
         for(i = 0; i < (n-1); i++)
         { 
            for(j = 0; j < n; j++){ S[i][j] = S[i+1][j];}
         }
         for(j = 0; j < n; j++){ S[n-1][j] = Scur[j];}
      }
      else 
      {
         for(j= 0; j < n; j++){ X[j] = Xold[j];}
         pGroup->WriteParams(X);
         m_pModel->SetObjFuncVal(Fcur);
         
         //restart directions
         for(i = 0; i < n; i++)
         { 
            for(j = 0; j < n; j++)
            { 
               S[i][j] = 0.00;
            }            
         }
         for(i = 0; i < n; i++) {S[i][i] = 1.00;}
         m_NumRestarts++;
      }
      m_AlgCount++;
      //end line minimization

      //write iteration result
      WriteRecord(m_pModel, iter+1, Fcur, fabs(Fold - Fcur));
      pStatus.pct = ((float)100.00*(float)(iter+1))/(float)m_MaxIter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //convergence termination criteria
      if(fabs(Fold - Fcur) <= Ftol)
      {
         count++;
         if(count >= MAX_COUNT)
         {
            pStatus.pct = 100.00;
            break;
         }
      }/* end if() */
      else
      {
         count = 0;
      }
      Fold = Fcur;

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */  

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write optimal results 
   WriteOptimal(m_pModel, Fcur);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //write algorithm metrics
   WriteAlgMetrics(this);

   //free up memory
   delete [] X;
   delete [] Xold;
   delete [] upr;
   delete [] lwr;  
   delete [] Scur;
   delete [] pmin;
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using Powell's algorithm.
******************************************************************************/
void PowellAlgorithm::Calibrate(void)
{ 
   FILE * pFile;
   char fileName[DEF_STR_SZ];
   int id;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);
     
   Optimize();

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   if(id == 0)
   {
      sprintf(fileName, "OstOutput%d.txt", id);

      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      m_pStats->WriteStats(stdout);
   }
} /* end Calibrate() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void PowellAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm        : Powell's Method (Conjugate Directions)\n");
   fprintf(pFile, "Max Iterations   : %d\n", m_MaxIter);
   fprintf(pFile, "Convergence Val  : %lf\n",m_ConvVal);
   fprintf(pFile, "Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Algorithm Evals  : %d\n", m_AlgCount);   
   //fprintf(pFile, "Total Evals      : %d\n", m_pModel->GetCounter());   
   fprintf(pFile, "Alg. Restarts    : %d\n", m_NumRestarts);
   fprintf(pFile, "Upper Violations : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations : %d\n", m_NumLwrViols);   
   m_pModel->WriteMetrics(pFile);
   m_pSearch->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
PWL_Program()

Calibrate or optimize using Powell's algorithm.
******************************************************************************/
void PWL_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("PowellAlgorithm", 1);
   PowellAlgorithm * PowAlg = new PowellAlgorithm(model);
   MEM_CHECK(PowAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ PowAlg->Calibrate(); }
   else { PowAlg->Optimize(); }

   delete PowAlg;
   delete model;
} /* end PWL_Program() */


