/******************************************************************************
File      : FletchReevesAlgorithm.cpp
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of the Fletcher-Reeves optimization algorithm.

Version History
08-26-03    lsm   created 
02-17-04    lsm   switched implementations
07-08-04    lsm   WriteSetup() is first action of algorithm
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
03-21-05    lsm   Added support for status file (used in grid computing)
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "FletchReevesAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "OptMathClass.h"
#include "OptSearchClass.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

#define EPS (NEARLY_ZERO)

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void FletchReevesAlgorithm::WarmStart(void)
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
FletchReevesAlgorithm::FletchReevesAlgorithm(ModelABC * pModel)
{
   FILE * inFile;
   ParameterGroup * pGroup;
   char * line;
   char tmp[DEF_STR_SZ];

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();   

   //init. everything to reasonable defaults
   m_MaxIter = 20;
   m_CurIter = 0;
   m_NumRestarts = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_AlgCount = 0;
   m_ConvVal = 1E-6; 
   m_MaxCount = MAX_COUNT;

   m_pStats = NULL;
   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   m_NumParams = pGroup->GetNumParams();

   NEW_PRINT("OptSearchClass", 1);
   m_pSearchAlg = new OptSearchClass(pModel);

   NEW_PRINT("OptMathClass", 1);
   m_pMath = new OptMathClass(pModel);
   
   MEM_CHECK(m_pMath);

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("FletchReevesAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginFletchReevesAlg", pFileName) == true)
   {
      FindToken(inFile, "EndFletchReevesAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginFletchReevesAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndFletchReevesAlg") == NULL)
      {
         if(strstr(line, "ConvergenceVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_ConvVal);
         }
         else if(strstr(line, "MaxStalls") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxCount);
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
void FletchReevesAlgorithm::Destroy(void)
{
   delete m_pSearchAlg;
   delete m_pMath;
   delete m_pStats;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Optimize()

Optimize the objective function using the Fletcher-Reeves algorithm.
******************************************************************************/
void FletchReevesAlgorithm::Optimize(void)
{
   StatusStruct pStatus;
   ParameterGroup * pGroup;
   double * p;
   int n, iter, count;
   double ftol, fret, old_fret;
   int j, its;
   double gg, gam, dgg, xmin, m, max, tmp, ftmp;
   double *g, *h, *xi, *upr, *lwr, * old_p, * pmin;
   const double * pGrad;

   n = m_NumParams;

   /* optimal parameter vector, can be updated by CalcGradient() */
   NEW_PRINT("double", n);
   pmin = new double[n];
   double fmin;

   NEW_PRINT("double", n);
   p = new double[n];

   NEW_PRINT("double", n);
   old_p = new double[n];

   NEW_PRINT("double", n);
   g = new double[n];

   NEW_PRINT("double", n);
   h = new double[n];

   NEW_PRINT("double", n);
   xi = new double[n];

   NEW_PRINT("double", n);
   upr = new double[n];

   NEW_PRINT("double", n);
   lwr = new double[n];
   
   MEM_CHECK(lwr);

   ftol = m_ConvVal;
   count = 0;

   //write setup to file
   WriteSetup(m_pModel, "Fletcher-Reeves");

   //initialize p
   pGroup = m_pModel->GetParamGroupPtr();
   pGroup->ReadParams(p);

   //handle warm start
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
      pGroup->ReadParams(p);
   }

   //get uppper and lower bounds
   for(j = 0; j < n; j++)
   {
      upr[j] = pGroup->GetParamPtr(j)->GetUprBnd();
      lwr[j] = pGroup->GetParamPtr(j)->GetLwrBnd();
   }/* end for() */

   //initialize fp and xi   
   fmin = old_fret = fret = m_pModel->Execute();
   pGroup->ReadParams(pmin);
   m_AlgCount++;

   //write banner and initial result to file
   WriteBanner(m_pModel, "iter  obj. function  ", "dObjFunc");
   WriteRecord(m_pModel, 0, fret, fret);
   //assemble status information and write to file
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxIter;
   pStatus.pct     = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //initialize the gradient
   its = 0;
   if(fret < fmin) fmin = fret;
   pGrad = m_pMath->CalcGradient(&fmin, pmin);
   if(fmin < fret) //found a better min during intermediate calculations
   {
      pGroup->WriteParams(pmin);
      pGroup->ReadParams(p);
      m_pModel->SetObjFuncVal(fmin);
      fret = fmin;
      old_fret = fret;
   }

   for(j = 0; j < n; j++){xi[j] = pGrad[j];}

   for(j = 0; j < n; j++)
   {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j];
   }
   
   for(its = its; its < m_MaxIter; its++)
   {
      if(IsQuit() == true){ break;}

      pStatus.curIter = m_CurIter = iter = its + 1;

      /*------------------------------------------------------
      normalize search direction so that largest value is 1.00
      -------------------------------------------------------*/
      //first find max value (in absolute terms)
      max = 0.00;
      for(j = 0; j < n; j++)
      { 
         tmp = fabs(xi[j]);
         if(tmp > max) { max = tmp;}
      }/* end for() */

      //now divide through by max
      for(j = 0; j < n; j++) { xi[j] = (xi[j] / max);}

      //line minimization
	  if(fret < fmin) fmin = fret;
      xmin = m_pSearchAlg->CalcStepSize(xi, &fmin, pmin);
	  
      for(j= 0; j < n; j++)
      {
         xi[j] *= xmin;
         old_p[j] = p[j];
         p[j] += xi[j];
         //if it's out of bounds, move half the distance to upr/lwr
         if(p[j] > upr[j]){p[j] = (upr[j]+old_p[j])/2.00; m_NumUprViols++;}
         if(p[j] < lwr[j]){p[j] = (old_p[j]+lwr[j])/2.00; m_NumLwrViols++;}
      }/* end for() */
      pGroup->WriteParams(p);
      ftmp = m_pModel->Execute();
      
      if(ftmp <= fret)
	  {
		  fret=ftmp;
	  }
      else //revert if moves didn't improve objective
      {
         pGroup->WriteParams(old_p);
		 pGroup->ReadParams(p);
         m_pModel->SetObjFuncVal(fret);
      }

	  if(fmin < fret) //intermediate step size found new optimal
	  {
		 pGroup->WriteParams(pmin);
	     pGroup->ReadParams(p);
		 m_pModel->SetObjFuncVal(fmin);
	     fret = fmin;
	  }

      m_AlgCount++;
      //end line minimization
      if(fret < fmin) fmin = fret;
      pGrad = m_pMath->CalcGradient(&fmin, pmin);
	  if(fmin < fret)
	  {
		 pGroup->WriteParams(pmin);
	     pGroup->ReadParams(p);
		 m_pModel->SetObjFuncVal(fmin);
	     fret = fmin;
	  }

      for(j = 0; j < n; j++){xi[j] = pGrad[j];}
      dgg=gg=0.00;

      for(j = 0; j < n; j++)
      {
         gg += g[j]*g[j];
         //dgg += xi[j]*xi[j]; Fletcher-Reeves
         dgg += (xi[j]+g[j])*xi[j]; //Polak-Ribiere
      }/* end for() */
      
      if(gg != 0.00)
      {
         gam = dgg/gg;
         for(j = 0; j < n; j++)
         {
            g[j] = -xi[j];
            xi[j] = h[j] = g[j] + gam*h[j];
         }/* end for() */
         //check slope   
         m = DotProduct(h, pGrad, n);      
         if(m >= 0) //use steepest descent
         {
            for(j = 0; j < n; j++)
            {
               xi[j] = h[j] = -pGrad[j];
               m_NumRestarts++;
            }
         }/* end if() */
      }/* end if() */

      //write iteration result to file
      WriteRecord(m_pModel, iter, fret, fabs(old_fret - fret));
      //assemble status information and write to file
      pStatus.pct  = ((float)100.00*(float)iter)/(float)m_MaxIter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //if gradient (gg) is 0.00, we are done
      if(gg == 0.00){pStatus.pct = 100.00; break;}

      //convergence termination criteria      
      if(fabs(old_fret - fret) <= ftol)
      {
         count++;
         if(count >= m_MaxCount)
         {
            pStatus.pct = 100.00; break;
         }
      }/* end if() */
      else
      {
         count = 0;
      }
      old_fret = fret;

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */      

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write results of final iteration to file
   WriteOptimal(m_pModel, fret);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   //write algorithm metrics
   WriteAlgMetrics(this);

   //free up memory
   delete [] p;
   delete [] g;
   delete [] h;
   delete [] xi;
   delete [] upr;
   delete [] lwr;
   delete [] old_p;
   delete [] pmin;
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using the Fletcher-Reeved algorithm.
******************************************************************************/
void FletchReevesAlgorithm::Calibrate(void)
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
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   if(id == 0)
   {
      sprintf(fileName, "OstOutput%d.txt", id);

      //write statistics of best parameter set to output file
      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      //write statistics of best parameter set to output file
      m_pStats->WriteStats(stdout);
   }/* end if() */   
} /* end Calibrate() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void FletchReevesAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm        : Fletcher-Reeves (Conjugate Gradient)\n");
   fprintf(pFile, "Max Iterations   : %d\n", m_MaxIter);
   fprintf(pFile, "Convergence Val  : %lf\n", m_ConvVal);
   fprintf(pFile, "Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Algorithm Evals  : %d\n", m_AlgCount);   
   //fprintf(pFile, "Total Evals      : %d\n", m_pModel->GetCounter());   
   fprintf(pFile, "S-D Restarts     : %d\n", m_NumRestarts);
   fprintf(pFile, "Upper Violations : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations : %d\n", m_NumLwrViols);   
   m_pModel->WriteMetrics(pFile);
   m_pMath->WriteMetrics(pFile);
   m_pSearchAlg->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
FLRV_Program()

Calibrate or optimize using the Fletcher-Reeves algorithm.
******************************************************************************/
void FLRV_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;
   
   NEW_PRINT("FletchReevesAlgorithm", 1);
   FletchReevesAlgorithm * MyAlg = new FletchReevesAlgorithm(model);
   
   MEM_CHECK(MyAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ MyAlg->Calibrate(); }
   else { MyAlg->Optimize(); }

   delete MyAlg;
   delete model;
} /* end FLRV_Program() */
