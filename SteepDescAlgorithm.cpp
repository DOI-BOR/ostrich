/******************************************************************************
File      : SteepDescAlgorithm.cpp
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of the Steepest Descent optimization algorithm.

Version History
08-26-03    lsm   created 
03-05-04    lsm   now reposrts obj. func. str.
07-08-04    lsm   WriteSetup() is now the first action of the algorithm
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
03-21-05    lsm   Added code to track algorithm status (used in grid computing)
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel. 
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SteepDescAlgorithm.h"
#include "Model.h"
#include "ModelBackup.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "OptSearchClass.h"
#include "OptMathClass.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"


/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void SteepDescAlgorithm::WarmStart(void)
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
SteepDescAlgorithm::SteepDescAlgorithm(ModelABC * pModel)
{
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   ParameterGroup * pGroup;

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();
   m_pStats = NULL;

   //init. everything to reasonable defaults
   m_MaxIter = 20;
   m_ConvVal = 1E-6; 
   m_AlgCount    = 0;
   m_CurIter     = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;

   
   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   m_NumParams = pGroup->GetNumParams();

   NEW_PRINT("double", m_NumParams);
   m_pSearchDir = new double[m_NumParams];   

   NEW_PRINT("OptSearchClass", 1);
   m_pSearchAlg = new OptSearchClass(pModel);   

   NEW_PRINT("OptMathClass", 1);
   m_pMath = new OptMathClass(pModel);
   MEM_CHECK(m_pMath);

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("SteepDescAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginSteepDescAlg", pFileName) == true)
   {
      FindToken(inFile, "EndSteepDescAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginSteepDescAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndSteepDescAlg") == NULL)
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
void SteepDescAlgorithm::Destroy(void)
{
   delete [] m_pSearchDir;

   delete m_pSearchAlg;
   delete m_pMath;
   delete m_pStats;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Optimize()

Optimize the objective function using the Steepest Descent algorithm.
******************************************************************************/
void SteepDescAlgorithm::Optimize(void)
{   
   StatusStruct pStatus;
   double curVal, oldVal, dObjFunc, step, tmp, max, upr, lwr, tst;
   int i, j;
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   Unchangeable1DArray pGrad;

   /* optimal parameter vector, can be updated by CalcGradient() */
   double * pmin;
   int n = m_NumParams;
   NEW_PRINT("double", n);
   pmin = new double[n];
   double fmin;

   //write setup to file
   WriteSetup(m_pModel, "Steepest Descent");

   //handle warm start
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }
   pGroup = m_pModel->GetParamGroupPtr();
   curVal = m_pModel->Execute();
   m_AlgCount++;
   dObjFunc = curVal;

   //write banner and initial result
   WriteBanner(m_pModel, "iter  obj. function  ", "dObjFunc");
   WriteRecord(m_pModel, 0, curVal, dObjFunc);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxIter;
   pStatus.pct = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //perform steepest descent iteration
   for(i = 0; i < m_MaxIter; i++)
   {
      if(IsQuit() == true){ break;}

      oldVal = curVal;
      pStatus.curIter = m_CurIter = i+1;

      //perform 1-D search using -gradient for search direction
     fmin = curVal;
     pGrad = m_pMath->CalcGradient(&fmin, pmin);
     //found a better min during gradient calculation?
     if(fmin < curVal)
     {
        for(j = 0; j < m_NumParams; j++)
        {
           pGroup->GetParamPtr(j)->SetEstVal(pmin[j]);
	    }
	    curVal = fmin;
	    m_pModel->SetObjFuncVal(fmin);
        oldVal = curVal;
     }

      for(j = 0; j < m_NumParams; j++){ m_pSearchDir[j] = -1.00 * pGrad[j];}

      /*------------------------------------------------------
      normalize search direction so that largest value is 1.00
      -------------------------------------------------------*/
      //first find max value (in absolute terms)
      max = 0.00;
      for(j = 0; j < m_NumParams; j++)
      { 
         tmp = fabs(m_pSearchDir[j]);
         if(tmp > max) { max = tmp;}
      }/* end for() */
      //now divide through by max, if max != 0
      if(fabs(max) > NEARLY_ZERO)
         for(j = 0; j < m_NumParams; j++) { m_pSearchDir[j] /= max; }

      //determine optimal step size and adjust search dir
      fmin = curVal;
      step = m_pSearchAlg->CalcStepSize(m_pSearchDir, &fmin, pmin);

      for(j = 0; j < m_NumParams; j++){ m_pSearchDir[j] *= step;}
      
      //make optimal move
      for(j = 0; j < m_NumParams; j++)
      {
         pParam = pGroup->GetParamPtr(j);
         tmp = pParam->GetEstVal();
         upr = pParam->GetUprBnd();
         lwr = pParam->GetLwrBnd();

         tst = tmp + m_pSearchDir[j];
         if(tst > upr){ tst = (tmp + upr)/2.00; m_NumUprViols++;}
         if(tst < lwr){ tst = (tmp + lwr)/2.00; m_NumLwrViols++;}
         pParam->SetEstVal(tst);
      }/* end for() */
      curVal = m_pModel->Execute();
      m_AlgCount++;

	  if(fmin < curVal) //found better minimum when computing step size
      {
         for(j = 0; j < m_NumParams; j++)
         {
            pGroup->GetParamPtr(j)->SetEstVal(pmin[j]);
		 }
		 m_pModel->SetObjFuncVal(fmin);
	     curVal = fmin;
	  }/* end if() */

      dObjFunc = fabs((oldVal - curVal)/(oldVal+NEARLY_ZERO));

      //write iteration result
      WriteRecord(m_pModel, (i+1), curVal, dObjFunc);
      pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_MaxIter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //converged?
      if(dObjFunc < m_ConvVal)
      {
         pStatus.pct = 100.00;
         break;
      }/* end if() */     

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false); 
   } /* end while() */

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write optimal result
   WriteOptimal(m_pModel, curVal);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //write algorithm metrics
   WriteAlgMetrics(this);

   delete [] pmin;
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using the Steepest Descent algorithm.
******************************************************************************/
void SteepDescAlgorithm::Calibrate(void)
{ 
   FILE * pFile;
   char fileName[DEF_STR_SZ];
   int id;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);
     
   Optimize();

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if(id == 0)
   {
      //compute statistics (variance and covariance)   
      m_pStats->CalcStats();

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
void SteepDescAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm        : Steepest-Descent\n");
   fprintf(pFile, "Max Iterations   : %d\n", m_MaxIter);
   fprintf(pFile, "Convergence Val  : %lf\n", m_ConvVal);
   fprintf(pFile, "Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Algorithm Evals  : %d\n", m_AlgCount);
   //fprintf(pFile, "Total Evals      : %d\n", m_pModel->GetCounter());
   fprintf(pFile, "Upper Violations : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations : %d\n", m_NumLwrViols);

   m_pModel->WriteMetrics(pFile);
   m_pMath->WriteMetrics(pFile);
   m_pSearchAlg->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
STPDSC_Program()

Calibrate or optimize using the steepest descent algorithm.
******************************************************************************/
void STPDSC_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("SteepDescAlgorithm", 1);
   SteepDescAlgorithm * SteepDesc = new SteepDescAlgorithm(model);
   MEM_CHECK(SteepDesc);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ SteepDesc->Calibrate(); }
   else { SteepDesc->Optimize(); }

   delete SteepDesc;
   delete model;
} /* end STPDSC_Program() */


