/******************************************************************************
File     : LevenbergAlgorithm.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The Levenberg-Marquardt algorithm is a hybrid numerical optimization method 
that initially uses the Steepest-Descent technique. However, since it is 
known that the Steepest-Descent algorithm converges very slowly near the 
optimum point, it is desirable to smoothly transition to a polynomial 
approximation method near the optimum. This implementation of the Levenberg 
algorithm is based upon the description/solution provided in the WinPEST 
user's manual, pages 9-42.

Version History
03-20-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-19-03    lsm   Parameter output uses WRITE_TX_BNR definition.
03-05-04    lsm   added obj. func. string to output, 
03-24-04    lsm   Added additional contstaints on parameter adjustments
07-08-04    lsm   WriteSetup() is first action of algorithm
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-11-05    lsm   Added algorithm metrics
03-21-05    lsm   Added support for status file (used in grid computing)
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics as well as 
                  the Jacobian can now be calculated in parallel. Added support
                  for the detection of temporarily insensitive parameters and
                  observations; when detected such insensitive terms are "held",
                  i.e. removed from consideration during the given iteration.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "LevenbergAlgorithm.h"
#include "ModelABC.h"
#include "Model.h"
#include "ModelBackup.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "StatsClass.h"

#include "Exception.h"
#include "Utility.h"
#include "WriteUtility.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void LevenbergAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
Calibrate()

Perform calibration using either Levenberg-Marquardt or GML-MS.
******************************************************************************/
void LevenbergAlgorithm::Calibrate(void)
{
   StatusStruct pStatus;
   char tmpStr[DEF_STR_SZ];
   char * gmlStr= (char *)"_GML";
   MyPoint Point, Best, Optimal;
   int i, j, id, np;
   double d, dmax, Fmin, init_lambda;
   m_NumEvals = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumMoveViols = 0;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   Point.v = NULL;
   Best.v = NULL;
   Optimal.v = NULL;

   //handle warm start
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }

   if(m_bMS == true)
   {
      init_lambda = m_Lambda;
      Fmin = NEARLY_HUGE;
      np = m_pModel->GetParamGroupPtr()->GetNumParams();
      Point.v = new double[np];
      Best.v = new double[np];
      Optimal.v = new double[np];

      //write setup
      WriteSetup(m_pModel, "GML-MS (multi-start #1)");
      
      for(i = 0; i < m_NumMS; i++)
      {
         SetIterationResidualsPrefix(gmlStr, 0); // set algorithm prefix
         SetIterationResidualsPrefix(NULL, i); // set trial prefix
         SetTrialNumber(i);

         m_Lambda = init_lambda; //reset lambda

         if(i != 0)
         {
            sprintf(tmpStr, "GML-MS (multi-start #%d)", i+1);
            WriteSetupNoDisclaimer(m_pModel, tmpStr);

            /* ------------------------------------------------------------
            Compute optimal new starting location> This is the point (out of
            several thousand trials) that is furthest from any previously 
            evaluated points.
            ------------------------------------------------------------ */
            dmax = 0.00;
            for(j = 0; j < 1000*np; j++)
            {
               //compute random starting location
               GetRndParamSet(&Point);
               d = GetMinDist(&Point, m_pList);
               if(d > dmax)
               {
                  dmax = d;
                  CopyPoint(&Point, &Best, np);
               }/* end if() */
            }/* end for() */
            //set new starting location
            m_pModel->GetParamGroupPtr()->WriteParams(Best.v);
         }/* end if() */
         
         CalibrateGML();

         if(m_Phi < Fmin)
         {
            Fmin = m_Phi;
            m_pModel->GetParamGroupPtr()->ReadParams(Optimal.v);
         }
      }/* end for() */
      //restore global optimal (best from all multi-starts)
      m_Phi = Fmin;
      m_pModel->GetParamGroupPtr()->WriteParams(Optimal.v);
   }/* end if() */
   else
   {
      //write setup
      WriteSetup(m_pModel, "Levenberg-Marquardt");

      CalibrateGML();
   }/* end else() */

   //write results of final iteration
   WriteOptimal(m_pModel, m_Phi);
   pStatus.numRuns = m_pModel->GetCounter();
   pStatus.curIter = m_NumIters;
   pStatus.maxIter = m_MaxIter;
   pStatus.pct = 100.00;
   WriteStatus(&pStatus);

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   {
      //write statistics to output file   
      char fileName[DEF_STR_SZ];
      FILE * pFile;

      sprintf(fileName, "OstOutput%d.txt", id);

      //write statistics of best parameter set to output file
      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      //write statistics of best parameter set to output file
      m_pStats->WriteStats(stdout);
   }/* end if() */

   //write algorithm metrics
   WriteAlgMetrics(this);

   delete [] Point.v;
   delete [] Best.v;
   delete [] Optimal.v;
}/* end Calibrate() */

/******************************************************************************
CopyPoint()

Copy data from one point to another.
******************************************************************************/
void LevenbergAlgorithm::CopyPoint(MyPoint * from, MyPoint * to, int np)
{
   int i;

   for(i = 0; i < np; i++)
   {
      to->v[i] = from->v[i];
   }
}/* end CopyPoint() */

/******************************************************************************
GetMinDist()

Compute minimum distance from prospective point to previously evaluated
points.
******************************************************************************/
double LevenbergAlgorithm::GetMinDist(MyPoint * Point, ParameterList * List)
{
   ParameterList * pCur;
   int i, np;
   double d, dmin;
   double v1, v2;

   dmin = NEARLY_HUGE;
   np = m_pModel->GetParamGroupPtr()->GetNumParams();

   for(pCur = List; pCur != NULL; pCur = pCur->pNxt)
   {
      d = 0.00;
      for(i = 0; i < np; i++)
      {
         v1 = Point->v[i];
         v2 = pCur->p.v[i];
         d += (v1-v2)*(v1-v2);
      }/* end for() */
      d = sqrt(d);

      if(d < dmin) dmin = d;
   }/* end for() */   

   return dmin;
}/* end GetMinDist() */

/******************************************************************************
GetRndParamSet()

Generate a random parameter set.
******************************************************************************/
void LevenbergAlgorithm::GetRndParamSet(MyPoint * Point)
{
   int i, np;
   double lwr, upr, range, r;

   np = m_pModel->GetParamGroupPtr()->GetNumParams();

   for(i = 0; i < np; i++)
   {
      lwr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLwrBnd();
      upr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUprBnd();
      range = upr - lwr;
      r = (double)MyRand() / (double)MY_RAND_MAX;
      Point->v[i] = (r * range) + lwr;
   }
}/* end GetRndParamSet() */

/******************************************************************************
CalibrateGML()

Perform calibration using Levenberg-Marquardt, equivalent to one application
of GML-MS.
******************************************************************************/
void LevenbergAlgorithm::CalibrateGML(void)
{
   static int GMLcount = 0;
   StatusStruct pStatus;
   double oldPhi;
   double * pMinJac;
   int done;
   int id;
   
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   m_CurIter = 0;

   //write banner
   WriteBanner(m_pModel, "iter  obj. function  ", "lambda");

   m_Phi = m_pModel->Execute();
   m_pModel->SaveBest(0); //save the input and output files of the best configuration
   m_BestSavedPhi = m_Phi;
   InsertParamSet();
   m_NumEvals++;
 
   //write iteration data
   WriteRecord(m_pModel, 0, m_Phi, m_Lambda);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxIter;
   pStatus.pct = (((float)100.00*(float)GMLcount)/(float)m_NumMS);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //main loop, iterates using Levenberg-Marquardt alg.
   done = 0;
   while(done == 0)
   {
      if(IsQuit() == true)
      { 
         done = 1;
         MPI_Bcast(&done, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
         break;
      }

      m_CurIter++;

      /* ----------------------------------------------------
      Calculate the Jacobian matrix, possibly in parallel.
      ------------------------------------------------------ */
      CalcJacobian();
      pMinJac = m_pStats->GetMinJac();

      if(id == 0)
      {
         CalcNormal(); //calculate the normal matrix
         CalcScale(); //calculate the scale matrix
         //determine best lambda for current iteration
         oldPhi = m_Phi;
         AdjustLambda();

         //check move against the best Jacobian evaluation
         if(pMinJac[0] < m_Phi)
         {
            m_Phi = pMinJac[0];
            //use the min Jacobian data to adjust model and parameter and observation groups
            m_pModel->SetObjFuncVal(pMinJac[0]);
            m_pModel->GetParamGroupPtr()->WriteParams(&(pMinJac[1]));
            m_pModel->GetObsGroupPtr()->WriteObservations(&(pMinJac[1+m_NumParams]));
         }/* end if() */
            
         //check for convergence
         m_PhiRatio = m_Phi / oldPhi;
         m_PhiRelRed = fabs(1.00 - m_PhiRatio);

         //write iteration data
         WriteRecord(m_pModel, m_CurIter, m_Phi, m_Lambda);
         pStatus.curIter = m_CurIter;
         pStatus.pct = (((float)100.00*(float)GMLcount)/(float)m_NumMS) + 
                       (((float)100.00*(float)m_CurIter)/((float)m_MaxIter*(float)m_NumMS));
         pStatus.numRuns = m_pModel->GetCounter();
         WriteStatus(&pStatus);

         if ((m_CurIter >= m_MaxIter) ||
            (m_PhiRelRed < m_Converge) ||
            ((oldPhi - m_Phi) < m_Converge))
         {
            done = 1;
            pStatus.pct = (((float)100.00*(float)GMLcount+1)/(float)m_NumMS);
         }/* end if() */
      }/* end if(master processor) */
      MPI_Bcast(&done, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end while() */
   m_NumIters = m_CurIter;
   GMLcount++;
}/* end CalibrateGML() */

/******************************************************************************
Optimize()

The Levenberg-Marquardt algorithm is a specialized algorithm that solves the
Least-Squares minimization problem. As such it is not suitable for general
optimization problems. The Optimize() routine has been coded as a no-op so 
that the inheritance from the Algorithm base class can be preserved.
******************************************************************************/
void LevenbergAlgorithm::Optimize(void)
{
   if(m_bMS == false)
   {
      printf("Levenberg-Marquardt algorithm can only be used for regression.\n");
   }
   else
   {
      printf("GML-MS algorithm can only be used for regression.\n");
   }
   return;
} /* end Optimize() */

/******************************************************************************
CalcJacobian()

Calculate the Jacobian matrix (matrix of dObs/dParams) using the algorithm
in the Statistics Class.
******************************************************************************/
void LevenbergAlgorithm::CalcJacobian(void)
{
   m_pJacob   = m_pStats->CalcJacobian(&m_BestSavedPhi);
   m_pJacobT  = m_pStats->GetJacobT();
   m_pJacobUW = m_pStats->GetJacobUW();
   m_pStats->AdjustJacobian();
}/* end CalcJacobian() */

/******************************************************************************
CalcNormal()

Calculate what is known as the 'normal' regression matrix:
   (J^T)*Q*J.
******************************************************************************/
void LevenbergAlgorithm::CalcNormal(void)
{   
   m_pNormal = m_pStats->CalcNormal();
}/* end CalcNormal() */

/******************************************************************************
CalcScale()

Calculate the scaling matrix. It is an all diagonal matrix that scales the
'normal' matrix so as to avoid numerical round-off errors and instability 
problems.
******************************************************************************/
void LevenbergAlgorithm::CalcScale(void)
{
   int i, p;   

   p = m_NumParams - m_pStats->GetNumHeldParams();
   for(i = 0; i < p; i++)
   {
      m_pScale[i][i] = (1.00 / sqrt(m_pNormal[i][i]));
   }/* end for() */  
}/* end CalcScale() */

/******************************************************************************
AdjustLambda()

Modify lambda in various ways to determine the best lambda for the current 
iteration. Reducing lambda is favored, but increasing lambda and leaving lambda
at its current value are also tested.
******************************************************************************/
void LevenbergAlgorithm::AdjustLambda(void)
{   
   int iter;
   double oldPhi, phiConst, phiDec, phiInc, phiTry;
   double lamConst, lamDec, lamInc, lamTry;   
      
   lamConst = m_Lambda;
   lamDec   = m_Lambda / m_LamSF;
   lamInc   = m_Lambda * m_LamSF;

   //display banner
   WriteInnerEval(WRITE_LEV, m_MaxLambdas, '.');
   
   //Compute initial lambda effects
   WriteInnerEval(1, m_MaxLambdas, '.');
   m_pInitBkup->Store();
   phiConst =  TryLambda(lamConst); //non-adjusted lambda trial
   m_pNonBkup->Store();

   WriteInnerEval(2, m_MaxLambdas, '-');
   m_pInitBkup->SemiRestore();
   phiDec = TryLambda(lamDec); //decreased lambda trial
   m_pDecBkup->Store(); 

   WriteInnerEval(3, m_MaxLambdas, '+');
   m_pInitBkup->SemiRestore();
   phiInc = TryLambda(lamInc); //increased lambda trial
   m_pIncBkup->Store();
   m_pInitBkup->SemiRestore();

   iter = 3;

   //check to see if none of the lambda adjustments were effective
   if((m_Phi < phiConst) && (m_Phi < phiDec) && (m_Phi < phiInc))
   {
      m_Lambda /= m_LamSF;
      WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'n');
      return;
   }/* end if() */

   /*------------------------------------------------
   Decreasing lambda caused obj. func. to decrease...
   and more so than a constant or increasing lambda
   ------------------------------------------------*/
   if((phiDec < m_Phi) && (phiDec <= phiInc) && (phiDec <= phiConst))
   {         
      oldPhi = m_Phi;
      m_pDecBkup->SemiRestore();
      while(iter < m_MaxLambdas)
      {
         //converged?
         m_PhiRatio = phiDec / oldPhi;
         m_PhiRelRed = 1.00 - m_PhiRatio;
         if((m_PhiRatio < m_RatioConv) || 
            ((m_PhiRelRed < m_RelRedConv) && (m_PhiRelRed > 0.00)))
         {
            //update phi, lambda
			m_pDecBkup->SemiRestore();
            m_Phi = phiDec;
            m_Lambda = lamDec;
            WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'c');
            return;
         }/* end if() */

         //try decreasing lambda
         WriteInnerEval(iter+1, m_MaxLambdas, '-');
         lamTry = lamDec / m_LamSF;
         oldPhi = phiDec;
         phiTry = TryLambda(lamTry);
         
         if(phiTry < oldPhi) //obj. reduced, accept trial move
         {   
            phiDec = phiTry;
            lamDec = lamTry;            
            m_pDecBkup->Store();            
         }/* end if() */
         else { break;} //failed to reduce, reject move and exit loop.
         iter++;         
      }/* end while */      
   }/* end if() */

   /*------------------------------------------------
   Increasing lambda caused obj. func. to decrease...
   and more so than a constant or decreasing lambda
   ------------------------------------------------*/
   if((phiInc < m_Phi) && (phiInc <= phiDec) && (phiInc <= phiConst))
   {         
      oldPhi = m_Phi;
      m_pIncBkup->SemiRestore();
      while(iter < m_MaxLambdas)
      {
         //converged?
         m_PhiRatio = phiInc / oldPhi;
         m_PhiRelRed = 1.00 - m_PhiRatio;
         if((m_PhiRatio < m_RatioConv) || 
            ((m_PhiRelRed < m_RelRedConv) && (m_PhiRelRed > 0.00)))
         {
            //update phi, lambda            
			m_pIncBkup->SemiRestore();
            m_Phi = phiInc;
            m_Lambda = lamInc;
            WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'c');
            return;
         }/* end if() */

         //try increasing lambda
         WriteInnerEval(iter+1, m_MaxLambdas, '+');
         lamTry = lamInc * m_LamSF;
         oldPhi = phiInc;
         phiTry = TryLambda(lamTry);
         
         if(phiTry < oldPhi) //obj. reduced, accept trial move
         {
            phiInc = phiTry;
            lamInc = lamTry;            
            m_pIncBkup->Store();            
         }/* end if() */
         else { break;}//failed to reduce, reject move and exit loop.
         iter++;         
      }/* end while */      
   }/* end if() */

   /*
   Didn't converge on a lambda, but some lambda(s) did reduce obj. func.
   Use the lambda that had the best result.
   */
   if((phiDec <= phiConst) && (phiDec <= phiInc))
   {     
      m_pDecBkup->SemiRestore();
      m_Phi = phiDec;
      m_Lambda = lamDec;
   }/* end if() */
   else if((phiConst <= phiDec) && (phiConst <= phiInc))
   {      
      m_pNonBkup->SemiRestore();
      m_Phi = phiConst;
      m_Lambda = lamConst;
   }/* end else if() */
   else
   {    
      m_pIncBkup->SemiRestore();
      m_Phi = phiInc;
      m_Lambda = lamInc;
   }/* end else() */

   WriteInnerEval(WRITE_ENDED, m_MaxLambdas, 'y');
}/* end AdjustLambda() */

/******************************************************************************
TryLambda()

Using the lambda argument, compute alpha (the Marquardt parameter) and the 
associated upgrade vector. Then apply this upgrade vector to the current set
of model parameters and execute the model.

Returns: The value of the obj. function, as returned by the Model class.
******************************************************************************/
double LevenbergAlgorithm::TryLambda(double lambda)
{
   double phi;

   //fill in vector of residuals [needed by CalcUpgrade() and CalcBeta()]
   m_pResid = m_pStats->CalcResiduals();
   m_pStats->AdjustResiduals();

   //solve for alpha (Marquardt parameter)
   CalcAlpha(lambda);
   //calculate the upgrade vector (the optimal direction)
   CalcUpgrade();
   //calculate the gamma vector (used in beta calculation)
   CalcGamma();
   /* Calculate beta, a factor that scales the upgrade vector 
   to achieve optimal step size */
   CalcBeta();
   //Update model parameters using beta and upgrade vector
   AdjModelParams();
   //Evaluate objective function at revised location
   phi = m_pModel->Execute();
   if(phi < m_BestSavedPhi)
   {
      m_pModel->SaveBest(0);
      m_BestSavedPhi = phi;
   }
   InsertParamSet();
   m_NumEvals++;
   return phi;
}/* end TryLambda() */

/******************************************************************************
CalcAlpha()

Compute alpha, the Marquardt parameter, based on the given lambda value. 
By definition, alpha = lambda/max(Si), where max(Si) is the maximum value of 
the diagonal elements of the scale matrix.
******************************************************************************/
void LevenbergAlgorithm::CalcAlpha(double lambda)
{
   int i, p;
   double max;
   double cur;

   max = m_pScale[0][0]*m_pScale[0][0];
   p = m_NumParams - m_pStats->GetNumHeldParams();
   for(i = 0; i < p; i++)
   {
      cur = m_pScale[i][i]*m_pScale[i][i];
      if(cur > max) {max = cur;}
   }/* end for() */
   
   m_Alpha = lambda/max;
}/* end CalcAlpha() */

/******************************************************************************
CalcUpgrade()

Compute the upgrade vector using a sequence of matrix and vector 
multiplication and inversion, as shown on page 20 of the WinPEST Manual.

Adjusted for Box-Cox Transformations ---- weights are incorporated into 
residuals prior to transformation.
******************************************************************************/
void LevenbergAlgorithm::CalcUpgrade(void)
{
   int i, n, p;

   p = m_NumParams - m_pStats->GetNumHeldParams();
   n = m_NumObs - m_pStats->GetNumHeldObs();

   //compute (JS)^T=S^T*J^T=S*J^T, where S^T = S, since S is a diag. matrix.
   //S is [pXp], J^T is [pXo], and S*J^T is [pXo]
   MatMult(m_pScale, m_pJacobT, m_pPbyO1, p, p, n);
   MatMult(m_pPbyO1, m_pJacob, m_pPbyP1, p, n, p);
   //multiply result by S [pxp], result is [pxp]
   MatMult(m_pPbyP1, m_pScale, m_pPbyP2, p, p, p);
   //add alpha*S^T*S to the result of previous step
   for(i = 0; i < p; i++)
   { m_pPbyP2[i][i] += (m_pScale[i][i]*m_pScale[i][i]*m_Alpha); }
   //invert the result of the previous step   
   MatInv(m_pPbyP2, m_pPbyP1, p);
   //multiply iverse by S*J^T, which is presently stored in m_pPbyO1
   MatMult(m_pPbyP1, m_pPbyO1, m_pPbyO2, p, p, n);
   VectMult(m_pPbyO2, m_pResid, m_pTmpVec, p, n);

   //premultiply result by S [pXp] and store in upgrade vector [pX1]
   VectMult(m_pScale, m_pTmpVec, m_pUpgrade, p, p);  
}/* end CalcUpgrade() */

/******************************************************************************
CalcGamma()

Gamma is used in the computation of beta, the optimum step size. Formulation 
of gamma is given on page 21 of the WinPEST Manual.
******************************************************************************/
void LevenbergAlgorithm::CalcGamma(void)
{
   int p, n;
   p = m_NumParams - m_pStats->GetNumHeldParams();
   n = m_NumObs - m_pStats->GetNumHeldObs();

   //gamma [ox1] = J [oxp] * u [px1]
   VectMult(m_pJacobUW, m_pUpgrade, m_pGamma, n, p);
}/* end CalcGamma() */

/******************************************************************************
CalcBeta()

Beta is the optimum step size for the direction specified by the upgrade 
vector. Formulation of beta is given on page 21 of the WinPEST Manual.
******************************************************************************/
void LevenbergAlgorithm::CalcBeta(void)
{
   int i, n;
   double numer;
   double denom;
   double wt;

   n = m_NumObs - m_pStats->GetNumHeldObs();
   numer = 0.00;
   denom = 0.00;

   for(i = 0; i < n; i++)
   {
      wt=GetObsWeight(m_pModel->GetObsGroupPtr()->GetObsPtr(i));
      numer += (m_pResid[i]*m_pGamma[i]*wt);
      denom += (m_pGamma[i]*m_pGamma[i]*wt*wt);
   }/* end for() */

   m_Beta = (numer/denom);   
}/* end CalcBeta() */

/******************************************************************************
AdjModelParams()

Modify model parameters using beta (step size) and upgrade vector (direction).
******************************************************************************/
void LevenbergAlgorithm::AdjModelParams(void)
{
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   double upr, lwr;
   double oldVal;
   double curVal;
   double adjust;
   double range;
   double maxAdj;
   int i;

   //must adjust the upgrade vector to account for held parameters
   m_pStats->AdjustVector(m_pUpgrade, false);

   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_NumParams; i++)
   {
      pParam = pGroup->GetParamPtr(i);
      oldVal = pParam->GetEstVal();
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      range = upr - lwr;
      maxAdj = range * m_MoveLimit;
      /*
      The optimal parameter adj. (not considering move limits) is 
      defined as: beta * (ugrade vector).
      */
      adjust = m_Beta * m_pUpgrade[i];

      /* 
      Check the optimal adjustment against move limits to prevent large 
      changes in parameters. We want to prevent large jumps because the 
      numerical solution relies on a Taylor Series expansion about the 
      current set of parameters. This is a linear approximation and 
      will only be valid in the proximity of the current paramter set.
      */      
      if(fabs(adjust) > maxAdj)
      {
         if(adjust > 0) {adjust = maxAdj;}
         else {adjust = -1.00 * maxAdj;}
         m_NumMoveViols++;
      } /* end if() */

      curVal = oldVal + adjust; 
   
      //if move exceeds limits, only go half the distance to the limit
      if(curVal <= lwr){ curVal = (oldVal+lwr)/2.00; m_NumLwrViols++;}
      if(curVal >= upr){ curVal = (oldVal+upr)/2.00; m_NumUprViols++;}

      pParam->SetEstVal(curVal);
   }/* end for() */   
}/* end AdjModelParams() */

/******************************************************************************
InsertParamSet()

Insert the most recently evaluated parameter set into the list.
******************************************************************************/
void LevenbergAlgorithm::InsertParamSet(void)
{
   if(m_bMS == false) return;

   ParameterList * pTmp;

   if(m_pList == NULL)
   {
      //allocate space for new list entry
      NEW_PRINT("ParameterList", 1);
      m_pList = new ParameterList;
      MEM_CHECK(m_pList);

      NEW_PRINT("double", m_NumParams);
      m_pList->p.v = new double[m_NumParams];
      MEM_CHECK(m_pList->p.v)

      m_pModel->GetParamGroupPtr()->ReadParams(m_pList->p.v);
      m_pList->pNxt = NULL;
   }/* end if(first entry) */
   else
   {
      for(pTmp = m_pList; pTmp->pNxt != NULL; pTmp = pTmp->pNxt)
      {
         //no-op, just advancing to end of list
      }

      //allocate space for new list entry
      NEW_PRINT("ParameterList", 1);
      pTmp->pNxt = new ParameterList;
      MEM_CHECK(pTmp->pNxt);
      pTmp = pTmp->pNxt;

      NEW_PRINT("double", m_NumParams);
      pTmp->p.v = new double[m_NumParams];
      MEM_CHECK(pTmp->p.v)

      m_pModel->GetParamGroupPtr()->ReadParams(pTmp->p.v);
      pTmp->pNxt = NULL;
   }/* end else(not first entry) */   
}/* end InsertParamSet() */

/******************************************************************************
CTOR()

Initialize data members to reasonable defaults. Certain defaults are then
overridden by the contents of the given confiugration fileName.
******************************************************************************/
LevenbergAlgorithm::LevenbergAlgorithm(ModelABC * pModel, bool bMulti)
{
   int i;
   int j;
   ParameterGroup * pParamGroup;
   ObservationGroup * pObsGroup;

   RegisterAlgPtr(this);

   m_Alpha  = 0.00;
   m_Beta   = 0.00;   
   m_Phi    = 1E6;
   m_PhiRatio = 1000.00;
   m_PhiRelRed = 1000.00;   
   m_Converge = 1E-4;
   m_RatioConv = 0.30;
   m_RelRedConv = 0.01;
   m_Lambda = 10.00;
   m_LamSF  = 1.1;
   m_MaxLambdas = 10;   
   m_MaxIter = 30;
   m_MoveLimit = 0.10; //limit parameter moves to 10% of current value
   m_bMS = bMulti;
   m_NumMS = 1;
   m_pList = NULL;

   m_pModel = pModel;
   pParamGroup = m_pModel->GetParamGroupPtr();
   pObsGroup = m_pModel->GetObsGroupPtr();

   m_NumObs = pObsGroup->GetNumObs();
   m_NumParams = pParamGroup->GetNumParams();

   NEW_PRINT("double", m_NumParams);
   m_pUpgrade = new double[m_NumParams];

   NEW_PRINT("double", m_NumParams);
   m_pTmpVec = new double[m_NumParams];

   NEW_PRINT("double", m_NumObs);
   m_pGamma = new double[m_NumObs];

   NEW_PRINT("double *", m_NumParams);
   m_pPbyP1 = new double *[m_NumParams];

   NEW_PRINT("double *", m_NumParams);
   m_pPbyP2 = new double *[m_NumParams];

   NEW_PRINT("double *", m_NumParams);
   m_pPbyO1 = new double*[m_NumParams];

   NEW_PRINT("double *", m_NumParams);
   m_pPbyO2 = new double*[m_NumParams];

   NEW_PRINT("double *", m_NumParams);
   m_pScale = new double*[m_NumParams];
   MEM_CHECK(m_pScale);

   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pPbyP1[i] = new double[m_NumParams];
   }
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pPbyP2[i] = new double[m_NumParams];
   }      
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumObs);
      m_pPbyO1[i] = new double[m_NumObs];
   }   
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumObs);
      m_pPbyO2[i] = new double[m_NumObs];
   }
   for(i = 0; i < m_NumParams; i++)
   {
      NEW_PRINT("double", m_NumParams);
      m_pScale[i] = new double[m_NumParams];
   }
   MEM_CHECK(m_pScale[i-1]);
   
   //fill in non-diagonal scale elements
   for(i = 0; i < m_NumParams; i++)
   {      
      for(j = 0; j < m_NumParams; j++)
      {         
         m_pScale[i][j] = 0.00;
      }/* end for() */
   }/* end for() */
            
   //setup model backups
   NEW_PRINT("ModelBackup", 1);
   m_pInitBkup = new ModelBackup(pModel);

   NEW_PRINT("ModelBackup", 1);
   m_pIncBkup = new ModelBackup(pModel);

   NEW_PRINT("ModelBackup", 1);
   m_pDecBkup = new ModelBackup(pModel);

   NEW_PRINT("ModelBackup", 1);
   m_pNonBkup = new ModelBackup(pModel);
   MEM_CHECK(m_pNonBkup);

   //set up statistics class   
   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);

   m_pJacob  = NULL;
   m_pJacobUW = NULL;
   m_pJacobT = NULL;
   m_pNormal = NULL;

   //configuration file can override certain defaults
   InitFromFile(GetInFileName());

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()

Frees up the various matrices and vectors used by the algorithm.
******************************************************************************/
void LevenbergAlgorithm::Destroy(void)
{
   ParameterList * pList, * pNext;
   int i;   

   delete [] m_pUpgrade;
   delete [] m_pTmpVec;
   delete [] m_pGamma;
 
   for(i = 0; i < m_NumParams; i++){delete [] m_pPbyP1[i];}
   delete [] m_pPbyP1;
   
   for(i = 0; i < m_NumParams; i++){delete [] m_pPbyP2[i];}
   delete [] m_pPbyP2;   
   
   for(i = 0; i < m_NumParams; i++){delete [] m_pPbyO1[i];}
   delete [] m_pPbyO1;
   
   for(i = 0; i < m_NumParams; i++){delete [] m_pPbyO2[i];}
   delete [] m_pPbyO2;

   for(i = 0; i < m_NumParams; i++){delete [] m_pScale[i];}
   delete [] m_pScale;

   //free up the model backups used in lambda trials
   delete m_pNonBkup;
   delete m_pInitBkup;
   delete m_pDecBkup;
   delete m_pIncBkup;
   delete m_pStats;

   if(m_pList != NULL)
   {
      pList = m_pList;
      while(pList->pNxt != NULL)
      {
         pNext = pList->pNxt; //isolate parameter
         pList->pNxt = pNext->pNxt; //unlink paramter
         //free up parameter
         delete [] pNext->p.v;
         delete pNext;
      }
      //free up head of list
      delete [] m_pList->p.v;
      delete m_pList;
   }

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void LevenbergAlgorithm::InitFromFile(IroncladString pLevFileName)
{   
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];

   //read in algorithm parameters
   pFile = fopen(pLevFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open algorithm config. file. Using Defaults");
      return;
   }/* end if() */   

   if(CheckToken(pFile, "BeginLevMar", pLevFileName) == true)
   {
      FindToken(pFile, "EndLevMar", pLevFileName);
      rewind(pFile);

      FindToken(pFile, "BeginLevMar", pLevFileName);
      line = GetNxtDataLine(pFile, pLevFileName);
      while(strstr(line, "EndLevMar") == NULL)
      {               
         //initial Marquardt lambda
         if(strstr(line, "InitialLambda") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Lambda); 
         }/*end if() */
         //Marquardt lambda scale factor
         else if(strstr(line, "LambdaScaleFactor") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_LamSF); 
         }/*end else if() */
         //parameter move limits
         else if(strstr(line, "MoveLimit") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_MoveLimit); 
         }/*end else if() */
         //convergence value
         else if(strstr(line, "AlgorithmConvergenceValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Converge); 
         }/*end else if() */
         //m_PhiRatio convergence value
         else if(strstr(line, "LambdaPhiRatio") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_RatioConv); 
         }/*end else if() */
         //m_PhiRelRed convergence value
         else if(strstr(line, "LambdaRelReduction") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_RelRedConv); 
         }/*end else if() */
         //maximum number of lambdas to try
         else if(strstr(line, "MaxLambdas") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxLambdas); 
         }/*end else if() */
         //maximum number of iterations
         else if(strstr(line, "MaxIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxIter); 
         }/*end else if() */
         //number of multi-starts
         else if(strstr(line, "NumMultiStarts") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumMS); 
            if(m_bMS == false) m_NumMS = 1;
            if(m_NumMS < 1) m_NumMS = 1;
         }/*end else if() */
         //intervals in zone of conv.
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */

         line = GetNxtDataLine(pFile, pLevFileName);
      } /* end while() */
   }/* end if() */   
   else
   {
      LogError(ERR_FILE_IO, "Using default algorithm setup.");
   }/* end else() */

   fclose(pFile);
}/* end InitFromFile() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void LevenbergAlgorithm::WriteMetrics(FILE * pFile)
{                      
   fprintf(pFile, "\nAlgorithm Metrics\n");
   if(m_bMS == false)
   {
      fprintf(pFile, "Algorithm         : Levenberg-Marquardt\n");
      fprintf(pFile, "Max Iterations    : %d\n", m_MaxIter);
      fprintf(pFile, "Actual Iterations : %d\n", m_NumIters);      
      fprintf(pFile, "Convergence Val   : %lf\n", m_Converge);
      fprintf(pFile, "LPRCV             : %lf\n", m_RatioConv);
      fprintf(pFile, "LRRCV             : %lf\n", m_RelRedConv);
      fprintf(pFile, "Max Lambda Trials : %d\n", m_MaxLambdas);
      fprintf(pFile, "Move Limit        : %lf\n", m_MoveLimit);
      fprintf(pFile, "Total Alg Evals   : %d\n", m_NumEvals);
      fprintf(pFile, "Total Evals       : %d\n", m_pModel->GetCounter());   
      fprintf(pFile, "Upper Violations  : %d\n", m_NumUprViols);
      fprintf(pFile, "Lower Violations  : %d\n", m_NumLwrViols);   
      fprintf(pFile, "Move Limit Viols  : %d\n", m_NumMoveViols); 
      fprintf(pFile, "LPRCV : Lambda-Phi Ratio Convergence Value\n");
      fprintf(pFile, "LRRCV : Lambda Relative Reduction Convergence Value\n");
   }
   else
   {
      fprintf(pFile, "Algorithm         : GML-MS\n");
      fprintf(pFile, "Num Multi-Starts  : %d\n", m_NumMS);
      fprintf(pFile, "Max Iterations    : %d (per multi-start)\n", m_MaxIter);
      fprintf(pFile, "Actual Iterations : %d (for last multi-start)\n", m_NumIters);
      fprintf(pFile, "Convergence Val   : %lf\n", m_Converge);
      fprintf(pFile, "LPRCV             : %lf\n", m_RatioConv);
      fprintf(pFile, "LRRCV             : %lf\n", m_RelRedConv);
      fprintf(pFile, "Max Lambda Trials : %d\n", m_MaxLambdas);
      fprintf(pFile, "Move Limit        : %lf\n", m_MoveLimit);
      fprintf(pFile, "Total Alg Evals   : %d (all multi-starts)\n", m_NumEvals);
      fprintf(pFile, "Avg Alg Evals     : %d (per multi-start)\n", m_NumEvals/m_NumMS);
      fprintf(pFile, "Total Evals       : %d (all multi-starts)\n", m_pModel->GetCounter());   
      fprintf(pFile, "Avg. Total Evals  : %d (per multi-start)\n", m_pModel->GetCounter()/m_NumMS);   
      fprintf(pFile, "Upper Violations  : %d (all multi-starts)\n", m_NumUprViols);
      fprintf(pFile, "Lower Violations  : %d (all multi-starts)\n", m_NumLwrViols);   
      fprintf(pFile, "Move Limit Viols  : %d (all multi-starts)\n", m_NumMoveViols); 
      fprintf(pFile, "Avg Upper Viols   : %d (per multi-start)\n", m_NumUprViols/m_NumMS);
      fprintf(pFile, "Avg Lower Viols   : %d (per multi-start)\n", m_NumLwrViols/m_NumMS);   
      fprintf(pFile, "Avg Mv Lmt Viols  : %d (per multi-start)\n", m_NumMoveViols/m_NumMS); 
      fprintf(pFile, "LPRCV : Lambda-Phi Ratio Convergence Value\n");
      fprintf(pFile, "LRRCV : Lambda Relative Reduction Convergence Value\n");
   }

   m_pModel->WriteMetrics(pFile);
   m_pStats->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
LEV_Program()

Create a model and solve using the Levenber-Marquardt algorithm.
******************************************************************************/
void LEV_Program(int argc, StringType argv[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;
   MEM_CHECK(model);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) 
   { 
      NEW_PRINT("LevenbergAlgorithm", 1);
      LevenbergAlgorithm * LA = new LevenbergAlgorithm(model, false);
      MEM_CHECK(LA);

      LA->Calibrate(); 
      delete LA;
   }
   else 
   { 
      printf("Levenberg-Marquardt algorithm can only be used for regression.\n");
   }/* end else() */
   
   delete model;
} /* end LEV_Program() */

/******************************************************************************
GMLMS_Program()

Create a model and solve using the Levenber-Marquardt algorithm with 
multi-starts.
******************************************************************************/
void GMLMS_Program(int argc, StringType argv[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;
   MEM_CHECK(model);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) 
   { 
      NEW_PRINT("LevenbergAlgorithm", 1);
      LevenbergAlgorithm * LA = new LevenbergAlgorithm(model, true);
      MEM_CHECK(LA);

      LA->Calibrate(); 
      delete LA;
   }
   else 
   { 
      printf("GML-MS algorithm can only be used for regression.\n");
   }/* end else() */
   
   delete model;
} /* end GMLMS_Program() */
