/******************************************************************************
ile      : SamplingAlgorithm.cpp
Author    : L. Shawn Matott
Copyright : 2007, L. Shawn Matott

An implementation of a sampling algorithm. Loosely based on Big Bang-Big Crunch
(BB-BC).

Version History
11-21-07    lsm   created 
******************************************************************************/
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "SamplingAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "LatinHypercube.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void SamplingAlgorithm::WarmStart(void)
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
SamplingAlgorithm::SamplingAlgorithm(ModelABC * pModel)
{
   ParameterGroup * pGroup;
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   int i;

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();

   //assign model and #params
   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   m_NumParams = pGroup->GetNumParams();   

   //init. everything to reasonable defaults
   m_bRndInit = true;
   m_MaxEvals = 100;
   m_MaxIter = 10;
   m_NumSamples = 10;
   m_Radius = 4.00;

   m_pAll = NULL;
   m_pSamples = NULL;
   m_pBest = NULL;
   m_pStats = NULL;
   m_pLwr = NULL;
   m_pUpr = NULL;
   m_pSD = NULL;
   m_pFwd = NULL;

   m_AlgCount = 0;
   m_CurIter = 0;

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("SamplingAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginSamplingAlg", pFileName) == true)
   {
      FindToken(inFile, "EndSamplingAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginSamplingAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndSamplingAlg") == NULL)
      {
         if(strstr(line, "MaxEvaluations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxEvals);
            if(m_MaxEvals < 1) m_MaxEvals = 1;
         }
         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default algorithm setup.");
   }/* end else() */

   fclose(inFile);

   m_MaxIter = (int)(sqrt((double)m_MaxEvals));
   m_NumSamples = m_MaxEvals/m_MaxIter;
   m_NumExtra = m_MaxEvals - (m_MaxIter*m_NumSamples);

   //allocate storage
   NEW_PRINT("MyPoint", (m_NumSamples+m_NumExtra));
   m_pSamples = new MyPoint[(m_NumSamples+m_NumExtra)];
   MEM_CHECK(m_pSamples);
   for(i = 0; i < (m_NumSamples+m_NumExtra); i++)
   {
      m_pSamples[i].F = NEARLY_HUGE;
      m_pSamples[i].ndim = m_NumParams;
      NEW_PRINT("double", m_NumParams);
      m_pSamples[i].v = new double[m_NumParams];
      MEM_CHECK(m_pSamples[i].v);
   }   

   NEW_PRINT("double", m_NumParams);
   m_pLwr = new double[m_NumParams];
   MEM_CHECK(m_pLwr);

   NEW_PRINT("double", m_NumParams);
   m_pUpr = new double[m_NumParams];
   MEM_CHECK(m_pUpr);
   
   for(i = 0; i < m_NumParams; i++)
   { 
      m_pLwr[i] = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLwrBnd();
      m_pUpr[i] = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUprBnd();
   }/* end for() */

   NEW_PRINT("double", m_NumParams);
   m_pSD = new double[m_NumParams];
   MEM_CHECK(m_pSD);
   for(i = 0; i < m_NumParams; i++)
   {      
      m_pSD[i] = 0.25*(m_pUpr[i] - m_pLwr[i]);      
   }

   NEW_PRINT("double", m_NumParams);
   m_pFwd = new double[m_NumParams];
   MEM_CHECK(m_pFwd);
   for(i = 0; i < m_NumParams; i++)
   {      
      m_pFwd[i] = 0.50;
   }

   IncCtorCount();
}/* end default CTOR */
      
/******************************************************************************
Destroy()

Frees up memory used by the algorithm.
******************************************************************************/
void SamplingAlgorithm::Destroy(void)
{
   int i;
   ParameterList * pList, * pNext;

   if(m_pAll != NULL)
   {
      pList = m_pAll;
      while(pList->pNxt != NULL)
      {
         pNext = pList->pNxt; //isolate parameter
         pList->pNxt = pNext->pNxt; //unlink paramter
         //free up parameter
         delete [] pNext->p.v;
         delete pNext;
      }
      //free up head of list
      delete [] m_pAll->p.v;
      delete m_pAll;
   }

   delete [] m_pSD;
   delete [] m_pLwr;
   delete [] m_pUpr;
   delete [] m_pFwd;

   if(m_pSamples != NULL)
   {
      for(i = 0; i < (m_NumSamples+m_NumExtra); i++)
      {
         delete [] m_pSamples[i].v;
      }
      delete [] m_pSamples;
   }
   delete m_pStats;   
   
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Optimize()

Optimize the objective function using the sampling algorithm.
******************************************************************************/
void SamplingAlgorithm::Optimize(void)
{   
   StatusStruct pStatus;
   int i, j, num_init;
   double F, F_old, dF;
   double r, sign, alpha, max_alpha;
   bool bWarmStart;

   m_CurIter = 0;

   //write setup
   WriteSetup(m_pModel, "Sampling Method");

   //write banner
   WriteBanner(m_pModel, "iter  obj. function  ", "relative change");

   //perform any user-requested evaluations
   UserDefinedEvaluations();
   
   // Determine the number of initial LHS samples
   bWarmStart = m_pModel->CheckWarmStart();
   if(bWarmStart == true)
   {
      num_init = 1;
   }
   else if(m_pModel->GetParamGroupPtr()->CheckExtraction() == true)
   {
      num_init = 1;
   }
   else if(m_bRndInit == true)
   {
      num_init = m_NumSamples+m_NumExtra;
   }
   else if(m_NumExtra > 0)
   {
      num_init = m_NumExtra;
      F_old = m_pBest->F;
   }
   else
   {
      num_init = 0;
   }

   // Generate and evaluate initial LHS samples
   if(num_init > 0)
   {
      GenerateInitialSamples(num_init);

      //handle warm start
      if(bWarmStart == true)
      {
         WarmStart();
         m_pModel->GetParamGroupPtr()->ReadParams(m_pSamples[0].v);
      }
      if(m_pModel->GetParamGroupPtr()->CheckExtraction() == true)
      {
         m_pModel->GetParamGroupPtr()->ReadParams(m_pSamples[0].v);
      }

       WriteInnerEval(WRITE_LHS, num_init, '.');
       for(i = 0; i < num_init; i++)
       {
          m_pModel->GetParamGroupPtr()->WriteParams(m_pSamples[i].v);   
          F = m_pModel->Execute();
          InsertParamSet(F);
          m_pSamples[i].F = F;
          m_AlgCount++;
          WriteInnerEval(i+1, 0, '.');
       }
       WriteInnerEval(WRITE_ENDED, 0, '.');

       /* -------------------------------------------------------------
       Report status
       -------------------------------------------------------------- */
       if(m_bRndInit == false)
       {
         if(m_pBest->F < F_old)
         {
            dF = fabs(100.00*((F_old - m_pBest->F)/F_old));
         }
         else
         {
            dF = 0.00;
         }
       }
       else
       {
          dF = 100.00;
       }       
       m_pModel->GetParamGroupPtr()->WriteParams(m_pBest->v); 
       WriteRecord(m_pModel, 0, m_pBest->F, dF);
       pStatus.curIter = 0;
       pStatus.maxIter = m_MaxIter;
       pStatus.pct = ((float)100.00*(float)0)/(float)m_MaxIter;
       pStatus.numRuns = m_pModel->GetCounter();
       WriteStatus(&pStatus);
   }/* end if() */
   
   /* --------------------------------------------------------
   Each iteration cooresponds to a random walk
   --------------------------------------------------------- */
   m_CurIter = 1;
   while(m_CurIter < m_MaxIter)
   {
      if(IsQuit() == true){ break;}
      
      //compute and evaluate random samples
      F_old = m_pBest->F;
      WriteInnerEval(WRITE_SMP, m_NumSamples, '.');
      for(i = 0; i < m_NumSamples; i++)
      {
         //generate random configuration
         for(j = 0; j < m_NumParams; j++)
         {
            //compute direction of displacement
            r = (double)MyRand()/(double)MY_RAND_MAX;
            if(r < m_pFwd[j])
            { 
               sign = +1.00;
               max_alpha = (m_pUpr[j] - m_pBest->v[j])/m_pSD[j]; // max positive displacement
            }
            else
            {
               sign = -1.00;
               max_alpha = (m_pBest->v[j] - m_pLwr[j])/m_pSD[j]; // max negative displacement
            }        
            alpha = m_Radius*((m_MaxIter - (m_CurIter - 1))/(double)m_MaxIter);
            if(alpha > max_alpha){ alpha = max_alpha;}
            r = (double)MyRand()/(double)MY_RAND_MAX;
            alpha *= r;
            m_pSamples[i].v[j] = m_pBest->v[j] + (sign*alpha*m_pSD[j]);
         }/*end for(each parameter) */

         m_pModel->GetParamGroupPtr()->WriteParams(m_pSamples[i].v);
         F = m_pModel->Execute();
         InsertParamSet(F);
         m_AlgCount++;
         WriteInnerEval(i+1, 0, '.');
      }/* end for(each sample) */
      WriteInnerEval(WRITE_ENDED, 0, '.');

      /* -------------------------------------------------------------
      Report status
      -------------------------------------------------------------- */
      dF = fabs(100.00*((F_old - m_pBest->F)/F_old));
      m_pModel->GetParamGroupPtr()->WriteParams(m_pBest->v); 
      WriteRecord(m_pModel, m_CurIter, m_pBest->F, dF);
      pStatus.curIter = m_CurIter;
      pStatus.maxIter = m_MaxIter;
      pStatus.pct = ((float)100.00*(float)m_CurIter)/(float)m_MaxIter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      m_CurIter++;

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for(outer iterations) */   

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write optimal results 
   m_pModel->GetParamGroupPtr()->WriteParams(m_pBest->v); 
   m_pModel->Execute();
   WriteOptimal(m_pModel, m_pBest->F);

   //write algorithm metrics
   WriteAlgMetrics(this);
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using the Sampling algorithm.
******************************************************************************/
void SamplingAlgorithm::Calibrate(void)
{ 
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);
     
   Optimize();

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   sprintf(fileName, "OstOutput0.txt");
   pFile = fopen(fileName, "a");   
   m_pStats->WriteStats(pFile);
   fclose(pFile);
   m_pStats->WriteStats(stdout);
} /* end Calibrate() */

/******************************************************************************
GenerateInitialSamples()

Create a set of LHS samples.
******************************************************************************/
void SamplingAlgorithm::GenerateInitialSamples(int num)
{
   int i, k;
   LatinHypercube * pLHS;
   pLHS = new LatinHypercube(m_NumParams, num);

   //prepare LHS sampler
   for(k = 0; k < m_NumParams; k++)
   { 
      pLHS->InitRow(k, m_pLwr[k], m_pUpr[k]);
   }/* end for() */

   //generate LHS samples
   for(i = 0; i < num; i++)
   {
      for(k = 0; k < m_NumParams; k++)
      {
         m_pSamples[i].v[k] = pLHS->SampleRow(k);
      }
   }

   //free up sampler
   delete pLHS;
}/* end GenerateInitialSamples() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void SamplingAlgorithm::WriteMetrics(FILE * pFile)
{
   int i;

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm            : Sampling Method\n");
   fprintf(pFile, "Max Evaluations      : %d\n", m_MaxEvals);
   fprintf(pFile, "Iterations           : %d\n", m_CurIter);
   fprintf(pFile, "Samples per Iter     : %d\n", m_NumSamples);
   fprintf(pFile, "User Defined Samples : %d\n", m_MaxEvals - (m_MaxIter*m_NumSamples) - m_NumExtra);
   fprintf(pFile, "Extra Samples        : %d\n", m_NumExtra);
   fprintf(pFile, "Algorithm Evals      : %d\n\n", m_AlgCount);

   fprintf(pFile, "\nParameter Standard Deviations (final estimate)\n");
   for(i = 0; i < m_NumParams; i++)
   {
      m_pModel->GetParamGroupPtr()->GetParamPtr(i)->Write(pFile, WRITE_BNR);
      fprintf(pFile, " : %E\n", m_pSD[i]);
   }/* end for() */

   fprintf(pFile, "\nParameter Forward Weights (at final iteration)\n");
   for(i = 0; i < m_NumParams; i++)
   {
      m_pModel->GetParamGroupPtr()->GetParamPtr(i)->Write(pFile, WRITE_BNR);
      fprintf(pFile, " : %E\n", m_pFwd[i]);
   }/* end for() */

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
UserDefinedEvaluations()

Evaluate objective functions of a user-defined list of parameter values.
Handy for seeding or restarting the SamplingAlgorithm from previous optimal.
******************************************************************************/
void SamplingAlgorithm::UserDefinedEvaluations(void)
{
   StatusStruct pStatus;
   ParameterList * pList;
   int j, k, count;
   char tmp[DEF_STR_SZ];
   const char * inFile = GetOstFileName();
   char * line, * pTok;
   FILE * pFile;

   pFile = fopen(inFile, "r");
   
   //skip routine if there aren't any user-defined configurations
   if(CheckToken(pFile, "BeginInitParams", inFile) == false)
   {
      fclose(pFile);
      return;
   }

   //verify proper tokens
   rewind(pFile);
   FindToken(pFile, "BeginInitParams", inFile);
   FindToken(pFile, "EndInitParams", inFile);
   rewind(pFile);

   //read in entries
   count = 0;
   FindToken(pFile, "BeginInitParams", inFile);
   line = GetNxtDataLine(pFile, inFile);
   while(strstr(line, "EndInitParams") == NULL)
   {
      pTok = line;

      //create new entry in list
      pList = InsertParamSet(NEARLY_HUGE);
      count++;

      //extract values, one-by-one, making any necessary conversions
      for(k = 0; k < m_NumParams; k++)
      {
         j = ExtractString(pTok, tmp);
         j = ValidateExtraction(j, k, m_NumParams, "SamplingAlgorithm::UserDefinedEvaluations()");
         pTok += j;
         pList->p.v[k] = m_pModel->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
      }/* end for() */                  

      line = GetNxtDataLine(pFile, inFile);
   }/* end while() */
   fclose(pFile);

   //perform model evaluations
   WriteInnerEval(WRITE_USR, count, '.');
   count = 0;
   for(pList = m_pAll; pList!= NULL; pList = pList->pNxt)
   {
      m_pModel->GetParamGroupPtr()->WriteParams(pList->p.v);   
      pList->p.F = m_pModel->Execute();
      m_AlgCount++;
      WriteInnerEval(++count, 0, '.');
   }
   WriteInnerEval(WRITE_ENDED, 0, '.');
   
   /* -------------------------------------------------
   Deduct user-defined evaluations from overall 
   sampling budget and revise configuration of sampler.
   ------------------------------------------------- */
   int max_evals;
   m_bRndInit = false;
   max_evals = m_MaxEvals - count;
   if(max_evals < 1) max_evals = 1;
   m_MaxIter = (int)(sqrt((double)max_evals));
   m_NumSamples = max_evals/m_MaxIter;
   m_NumExtra = max_evals - (m_MaxIter*m_NumSamples);
   
   /* --------------------------------------------------------------
   Report status
   --------------------------------------------------------------- */
   m_pModel->GetParamGroupPtr()->WriteParams(m_pBest->v);
   WriteRecord(m_pModel, -1, m_pBest->F, 100.00);
   pStatus.curIter = -1;
   pStatus.maxIter = m_MaxIter;
   pStatus.pct = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
}/* end UserDefinedEvaluations() */

/******************************************************************************
InsertParamSet()

Insert the most recently evaluated parameter set into the list, preserving 
order of list.
******************************************************************************/
ParameterList * SamplingAlgorithm::InsertParamSet(double F)
{
   int j;
   ParameterList * pTmp, * pPrev;

   if(m_pAll == NULL)
   {
      //allocate space for new list entry
      NEW_PRINT("ParameterList", 1);
      m_pAll = new ParameterList;
      MEM_CHECK(m_pAll);

      NEW_PRINT("double", m_NumParams);
      m_pAll->p.v = new double[m_NumParams];
      MEM_CHECK(m_pAll->p.v)

      m_pModel->GetParamGroupPtr()->ReadParams(m_pAll->p.v);
      m_pAll->p.F = F;
      m_pAll->pNxt = NULL;

      m_pBest = &(m_pAll->p);

      for(j = 0; j < m_NumParams; j++)
      {
         m_pSD[j] = m_pBest->v[j]*((m_MaxIter - m_CurIter)/(double)m_MaxIter);
      }/* end for() */

      return m_pAll;
   }/* end if(first entry) */
   else
   {
      //locate correct position in list to preserve increasing order of F
      pPrev = NULL;
      for(pTmp = m_pAll; pTmp != NULL; pTmp = pTmp->pNxt)
      {
         if(F < pTmp->p.F)
         {
            break;
         }
         pPrev = pTmp;
      }/* end for() */

      if(pPrev == NULL) //insert at head of list
      {
         //allocate space for new list entry
         NEW_PRINT("ParameterList", 1);
         m_pAll = new ParameterList;
         MEM_CHECK(m_pAll);

         NEW_PRINT("double", m_NumParams);
         m_pAll->p.v = new double[m_NumParams];
         MEM_CHECK(m_pAll->p.v)

         m_pModel->GetParamGroupPtr()->ReadParams(m_pAll->p.v);
         m_pAll->p.F = F;
         m_pAll->pNxt = pTmp;

         m_pBest = &(m_pAll->p);

         for(j = 0; j < m_NumParams; j++)
         {
            m_pSD[j] = m_pBest->v[j]*((m_MaxIter - m_CurIter)/(double)m_MaxIter);

            if(m_pBest->v[j] < pTmp->p.v[j]) //reducing parameter improved obj
            {
               m_pFwd[j] = 0.00;
            }
            else //increasing parameter improved obj
            {
               m_pFwd[j] = 1.00;
            }
         }/* end for() */
         return m_pAll;
      }
      else
      {
         //allocate space for new list entry
         NEW_PRINT("ParameterList", 1);
         pPrev->pNxt = new ParameterList;
         MEM_CHECK(pPrev->pNxt);
         pPrev = pPrev->pNxt;

         NEW_PRINT("double", m_NumParams);
         pPrev->p.v = new double[m_NumParams];
         MEM_CHECK(pPrev->p.v)

         m_pModel->GetParamGroupPtr()->ReadParams(pPrev->p.v);
         pPrev->p.F = F;
         pPrev->pNxt = pTmp;

         //reverse direction
         for(j = 0; j < m_NumParams; j++)
         {
            m_pFwd[j] = 1.00 - m_pFwd[j];
         }
         return pPrev;
      }/* end else(not first entry) */   
   }/* end if() */
}/* end InsertParamSet() */

/******************************************************************************
SMP_Program()

Calibrate or optimize using sampling algorithm.
******************************************************************************/
void SMP_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("SamplingAlgorithm", 1);
   SamplingAlgorithm * SmpAlg = new SamplingAlgorithm(model);
   MEM_CHECK(SmpAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ SmpAlg->Calibrate(); }
   else { SmpAlg->Optimize(); }

   delete SmpAlg;
   delete model;
} /* end SMP_Program() */


