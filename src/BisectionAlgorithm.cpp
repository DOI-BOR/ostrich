/******************************************************************************
ile      : BisectionAlgorithm.cpp
Author    : L. Shawn Matott
Copyright : 2007, L. Shawn Matott

An implementation of a basic bisection algorithm.

Version History
08-01-07    lsm   created 
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "BisectionAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void BisectionAlgorithm::WarmStart(void)
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
BisectionAlgorithm::BisectionAlgorithm(ModelABC * pModel)
{
   ParameterGroup * pGroup;
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();
   m_pStats = NULL;

   //init. everything to reasonable defaults
   m_MaxOuter = 50;
   m_MaxInner = 20;
   m_ConvVal = 1E-6; 
   m_AlgCount = 0;

   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   m_NumParams = pGroup->GetNumParams();   

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("BisectionAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginBisectionAlg", pFileName) == true)
   {
      FindToken(inFile, "EndBisectionAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginBisectionAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndBisectionAlg") == NULL)
      {
         if(strstr(line, "ConvergenceVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_ConvVal);
         }
         else if(strstr(line, "MaxOuterIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxOuter);
         }
         else if(strstr(line, "MaxInnerIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxInner);
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
void BisectionAlgorithm::Destroy(void)
{
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Optimize()

Optimize the objective function using the bisection algorithm.
******************************************************************************/
void BisectionAlgorithm::Optimize(void)
{
   double Xcur, Xupr, Xlwr, Xmid, Xqtr, X3qt, Xmin;
   double Fcur, Fupr, Flwr, Fmid, Fqtr, F3qt, Fmin, Fold;
   double * pXmin;
   double range, r, val;

   StatusStruct pStatus;
   int i, j, p, tmp;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   tmp = 0;

   //write setup
   WriteSetup(m_pModel, "Bisection Method");

   //read in best result from previous run, if desired
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }

   pXmin = new double[m_NumParams];
   Fold = Fcur = m_pModel->Execute();
   m_pModel->GetParamGroupPtr()->ReadParams(pXmin);
   m_AlgCount++;

   //write banner and initial result
   WriteBanner(m_pModel, "iter  obj. function  ", "dObjFunc");
   WriteRecord(m_pModel, 0, Fcur, Fcur);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxOuter;
   pStatus.pct = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   /* ---------------------------------
   One outer iteration cooresponds to
   application bisection to each 
   parameter.
   --------------------------------- */
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_MaxOuter; i++)
   {
      if(IsQuit() == true){ break;}

      if(i > 0)
      {
         for(j = 0; j < m_NumParams; j++) //for each parameter
         {
            //generate a random between lower and upper bound
            Xlwr = pGroup->GetParamPtr(j)->GetLwrBnd();
            Xupr = pGroup->GetParamPtr(j)->GetUprBnd();
            range = Xupr - Xlwr;
            r = (double)MyRand() / (double)MY_RAND_MAX;
            val = (r * range) + Xlwr;
            pGroup->GetParamPtr(j)->SetEstVal(val);
         }/* end for() */
      }/* end if() */

      pStatus.curIter = m_CurIter = i+1;

      WriteInnerEval(WRITE_BIS, 5*m_NumParams+m_MaxInner*m_NumParams*2, '.');
      tmp =0;
      for(p = 0; p < m_NumParams; p++)
      {
         pParam = pGroup->GetParamPtr(p);
         Xcur = pParam->GetEstVal();

         /*------------------------------------------
         Assign and evaluate initial points, these
         subdivide the domain into 4 quadrants.
         ----------------------------------------- */
         Xupr = pParam->GetUprBnd();
         pParam->SetEstVal(Xupr);
         Fupr = m_pModel->Execute();
         m_AlgCount++;
         WriteInnerEval(++tmp, 0, '.');

         Xlwr = pParam->GetLwrBnd();
         pParam->SetEstVal(Xlwr);
         Flwr = m_pModel->Execute();
         m_AlgCount++;
         WriteInnerEval(++tmp, 0, '.');

         Xqtr = Xlwr + 0.25*(Xupr - Xlwr);
         pParam->SetEstVal(Xqtr);
         Fqtr = m_pModel->Execute();
         m_AlgCount++;
         WriteInnerEval(++tmp, 0, '.');

         Xmid = Xlwr + 0.50*(Xupr - Xlwr);
         pParam->SetEstVal(Xmid);
         Fmid = m_pModel->Execute();
         m_AlgCount++;
         WriteInnerEval(++tmp, 0, '.');

         X3qt = Xlwr + 0.75*(Xupr - Xlwr);
         pParam->SetEstVal(X3qt);
         F3qt = m_pModel->Execute();
         m_AlgCount++;
         WriteInnerEval(++tmp, 0, '.');

         /* ---------------------------------------------
         Perform bisections. Each bisection will reduce
         the search space by 50%.
         --------------------------------------------- */
         for(j = 0; j < m_MaxInner; j++)
         {            
            //is mid-point current minimum?
            if((Fmid <= Fupr) && (Fmid <= Flwr) && (Fmid <= Fqtr) && (Fmid <= F3qt))
            {
               //store new minimum
               Xmin = Xmid;
               Fmin = Fmid;

               //shrink domain
               Xlwr = Xqtr;
               Flwr = Fqtr;
               Xupr = X3qt;
               Fupr = F3qt;
            }/* end if() */
            //is quarter-point current minimum?
            else if((Fqtr <= Fupr) && (Fqtr <= Flwr) && (Fqtr <= Fmid) && (Fqtr <= F3qt))
            {
               //store new minimum
               Xmin = Xqtr;
               Fmin = Fqtr;         

               //shrink domain
               Xupr = Xmid;
               Fupr = Fmid;
               Xmid = Xqtr;
               Fmid = Fqtr;         
            }/* end if() */
            //is three-quarter-point current minimum?
            else if((F3qt <= Fupr) && (F3qt <= Flwr) && (F3qt <= Fmid) && (F3qt <= Fqtr))
            {
               //store new minimum
               Xmin = X3qt;
               Fmin = F3qt;
      
               //shrink domain
               Xlwr = Xmid;
               Flwr = Fmid;
               Xmid = X3qt;
               Fmid = F3qt;
            }/* end if() */
            //is upper-bound current minimum?
            else if((Fupr <= F3qt) && (Fupr <= Flwr) && (Fupr <= Fmid) && (Fupr <= Fqtr))
            {
               //store the new minimum
               Xmin = Xupr;
               Fmin = Fupr;
         
               //shrink domain
               Xlwr = X3qt;
               Flwr = F3qt;

               //assign and evaluate new mid-point
               Xmid = Xlwr + 0.5*(Xupr - Xlwr);
               pParam->SetEstVal(Xmid);
               Fmid = m_pModel->Execute();
               m_AlgCount++;
               WriteInnerEval(++tmp, 0, '.');
            }/* end if() */
            //is lower-bound current minimum?
            else if((Flwr <= F3qt) && (Flwr <= Fupr) && (Flwr <= Fmid) && (Flwr <= Fqtr))
            {
               //store the new minimum
               Xmin = Xlwr;
               Fmin = Flwr;

               //shrink domain
               Xupr = Xqtr;
               Fupr = Fqtr;

               //assign and evaluate new mid-point
               Xmid = Xlwr + 0.5*(Xupr - Xlwr);
               pParam->SetEstVal(Xmid);
               Fmid = m_pModel->Execute();
               m_AlgCount++;
               WriteInnerEval(++tmp, 0, '.');
            }/* end if() */
            //assume mid-point
            else
            {
               //store new minimum
               Xmin = Xmid;
               Fmin = Fmid;
               //fprintf(pFile, "unk ");

               //shrink domain
               Xlwr = Xqtr;
               Flwr = Fqtr;
               Xupr = X3qt;
               Fupr = F3qt;
            }/* end if() */

            //assign and evaluate new quarter points
            Xqtr = Xlwr + 0.25*(Xupr - Xlwr);
            pParam->SetEstVal(Xqtr);
            Fqtr = m_pModel->Execute();
            m_AlgCount++;
            WriteInnerEval(++tmp, 0, '.');

            X3qt = Xlwr + 0.75*(Xupr - Xlwr);
            pParam->SetEstVal(X3qt);
            F3qt = m_pModel->Execute();
            m_AlgCount++;
            WriteInnerEval(++tmp, 0, '.');
         }/* end for(inner iterations) */
         if(Fmin < Fcur)
         {
            Fcur = Fmin;
            pParam->SetEstVal(Xmin);
            m_pModel->SetObjFuncVal(Fmin);
            m_pModel->GetParamGroupPtr()->ReadParams(pXmin);
         }
         else //search failed revert to initial value
         {
            pParam->SetEstVal(Xcur);
            m_pModel->SetObjFuncVal(Fcur);
         }
      }/* end for(parameters) */
      WriteInnerEval(WRITE_ENDED, 0, '.');

      if(Fold < Fmin)
      {
         m_pModel->SetObjFuncVal(Fold);
         m_pModel->GetParamGroupPtr()->WriteParams(pXmin);
      }

      //write iteration result
      WriteRecord(m_pModel, i+1, Fcur, fabs(Fold - Fcur));
      pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_MaxOuter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      
      Fold = Fcur;

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for(outer iterations) */   

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write optimal results 
   WriteOptimal(m_pModel, Fcur);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //write algorithm metrics
   WriteAlgMetrics(this);

   delete [] pXmin;
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using Bisection algorithm.
******************************************************************************/
void BisectionAlgorithm::Calibrate(void)
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
void BisectionAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm        : Bisection Method\n");
   fprintf(pFile, "Max Outer Iters  : %d\n", m_MaxOuter);
   fprintf(pFile, "Max Inner Iters  : %d\n", m_MaxInner);
   fprintf(pFile, "Convergence Val  : %lf\n",m_ConvVal);
   fprintf(pFile, "Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Algorithm Evals  : %d\n", m_AlgCount);   
   //fprintf(pFile, "Total Evals      : %d\n", m_pModel->GetCounter());   
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
BIS_Program()

Calibrate or optimize using Bisection algorithm.
******************************************************************************/
void BIS_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("BisectionAlgorithm", 1);
   BisectionAlgorithm * BisAlg = new BisectionAlgorithm(model);
   MEM_CHECK(BisAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ BisAlg->Calibrate(); }
   else { BisAlg->Optimize(); }

   delete BisAlg;
   delete model;
} /* end BIS_Program() */


