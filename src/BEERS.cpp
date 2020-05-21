/******************************************************************************
File     : BEERS.cpp
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

Balanced Exploration-Exploitation Random Search (BEERS) applies a balanced
approach to randomly search a parameter space. The search initially favors
exploration of the search space and will try to maximize the distance between 
points that are evaluated. As the search progresses the algorithm favors 
exploitation and will favor points that are close to the current optimal. To
facilitate exploration vs. exploitation the algorithm maintains an archive of
every point that is ever evaluated during the search.

Version History
05-16-15    lsm   added copyright information and initial comments.
******************************************************************************/

#include "mpi_stub.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "BEERS.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "Utility.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void BEERS::WarmStart(void)
{
   int i;
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);

   for(i = 0; i < np; i++)
   {
      m_pArchive->X[i] = pbest[i];
   }
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
BEERS::BEERS(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pArchive = NULL;
   m_pBest = NULL;
   m_pMin = NULL;
   m_pMax = NULL;
   m_pRange = NULL;
   m_NumSamples = 0;
   m_CurSample = 0;
   m_MinProbAccept = 0;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the PSO and it's member variables.
******************************************************************************/
void BEERS::Destroy(void)
{
   m_NumSamples = 0;
   m_CurSample = 0;
   delete [] m_pMin;
   delete [] m_pMax;
   delete [] m_pRange;
   DestroyArchive(m_pArchive);
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using BEERS.
******************************************************************************/
void BEERS::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using BEERS.
******************************************************************************/
void BEERS::Optimize(void)
{
   StatusStruct pStatus;
   ParameterGroup * pGroup;
   int num;
   bool bBanner = false;

   InitFromFile(GetInFileName());

   WriteSetup(m_pModel, "BEERS - Balanced Exploration-Exploitation Random Search");
   //write banner
   WriteBanner(m_pModel, "iter   best_value     ", "Samples_Remaining");

   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();

   m_pMin = new double[num];
   m_pMax = new double[num];
   m_pRange = new double[num];

   NEW_PRINT("ArchiveStruct", 1);
   m_pArchive = new ArchiveStruct;
   MEM_CHECK(m_pArchive);

   m_pArchive->nF = 1;
   m_pArchive->nX = num;
   m_pArchive->F = new double[1];
   m_pArchive->F[0] = NEARLY_HUGE;
   m_pArchive->X = new double[num];
   m_pArchive->pNext = NULL;
   m_pArchive->Z = 0; //rank

   double lwr, upr, range, r, rval;
   for(int i = 0; i < num; i++)
   {
      //generate a random between lower and upper bound
      lwr = pGroup->GetParamPtr(i)->GetLwrBnd();
      upr = pGroup->GetParamPtr(i)->GetUprBnd();
      range = upr - lwr;
      r = (double)MyRand() / (double)MY_RAND_MAX;
      rval = (r * range) + lwr;
      m_pArchive->X[i] = rval;

      m_pMin[i] = lwr;
      m_pMax[i] = upr;
      m_pRange[i] = range;
   }/* end for() */

   //read in best result from previous run, if desired
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }
   // extract values from existing output file, if desired
   if(pGroup->CheckExtraction() == true)
   {
      pGroup->ReadParams(m_pArchive->X);
   } 

   //evaluate first sample
   WriteInnerEval(WRITE_SMP, m_NumSamples, '.');
   WriteInnerEval(1, m_NumSamples, '.');
   pGroup->WriteParams(m_pArchive->X);
   m_pArchive->F[0] = m_pModel->Execute();
   m_pBest = m_pArchive;
   WriteInnerEval(WRITE_ENDED, m_NumSamples, '.');

   //write initial config.
   WriteRecord(m_pModel, 0, m_pBest->F[0], m_NumSamples-1);
   m_CurSample = pStatus.curIter = 1;
   pStatus.maxIter = m_NumSamples;
   pStatus.pct = (float)100/(float)m_NumSamples;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //main optimization loop   
   ArchiveStruct * pCur;
   ArchiveStruct * pTail = m_pArchive;
   bool bCandidateAccepted;
   double Pexploit, Pexplore;
   double Wexploit, Wexplore;
   double Paccept, Raccept;
   WriteInnerEval(WRITE_SMP, m_NumSamples, '.');
   m_MaxDist = sqrt((double)num);
   for(int g = 1; g < m_NumSamples; g++)
   {      
      if(IsQuit() == true){ break;}

      //estimate the maximum distance between a candidate point and the nearest point in the archive
      //m_MaxDist = EstimateMaxDistance(m_pArchive, m_pMin, m_pRange);

      //Assign model probabilities
      AssignModelProbs(m_pArchive, m_pBest->F[0]);

      WriteArchive();

      //compute exploration and exploitation weights
      Wexploit = ((double)(g+1) / (double)m_NumSamples);
      Wexplore = 1.00 - Wexploit;

      //compute new candidate solution
      pCur = new ArchiveStruct;
      pCur->nF = 1;
      pCur->nX = num;
      pCur->F = new double[1];
      pCur->F[0] = NEARLY_HUGE;
      pCur->X = new double[num];
      pCur->pNext = NULL;
      pCur->Z = 0;

      bCandidateAccepted = false;
      m_MinProbAccept = 0.00;
      while(bCandidateAccepted == false)
      {
         //generate random candidate
         for(int i = 0; i < num; i++)
         {
            pCur->X[i] = (((double)MyRand() / (double)MY_RAND_MAX) * m_pRange[i]) + m_pMin[i];
         }/* end for() */

         //compute exploitation and exploration probabilities
         CalcProbabilities(pCur, m_pArchive, m_MaxDist, m_pMin, m_pRange, &Pexploit, &Pexplore);

         Paccept = (Wexploit * Pexploit) + (Wexplore * Pexplore);
         if(Paccept < m_MinProbAccept)
         {
            Paccept = m_MinProbAccept;
         }

         Raccept = ((double)MyRand() / (double)MY_RAND_MAX);
         if(Raccept < Paccept)
         {
            bCandidateAccepted = true;
         }/* end if() */

         //ensures that a candidate is eventually selected
         m_MinProbAccept += 1E-6;
      }/* end while() */
      
      //evaluate parameter set
      pGroup->WriteParams(pCur->X);
      pCur->F[0] = m_pModel->Execute();
      m_CurSample = pStatus.curIter = g+1;
      WriteInnerEval(g+1, m_NumSamples, '.');
      bBanner = false;

      //add to archive
      pTail->pNext = pCur;
      pTail = pCur;

      //update best solution
      if(pCur->F[0] < m_pBest->F[0])
      {
         m_pBest = pCur;
         WriteInnerEval(WRITE_ENDED, m_NumSamples, '.');
         WriteRecord(m_pModel, (g+1), m_pBest->F[0], m_NumSamples-g-1);
         if(g < m_NumSamples)
            WriteInnerEval(WRITE_SMP, m_NumSamples, '.');
         bBanner = true;                  
      }/* end if() */

      AdjustRanks(pCur, m_pArchive);

      pStatus.pct = ((float)100.00*(float)(g+1))/(float)m_NumSamples;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
   }/* end for() */

   if(bBanner ==  false)
      WriteInnerEval(WRITE_ENDED, 0, '.');

   pGroup->WriteParams(m_pBest->X);
   WriteOptimal(m_pModel, m_pBest->F[0]);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //write algorithm metrics
   WriteAlgMetrics(this);
} /* end Optimize() */

/******************************************************************************
AdjustRanks()

Iterate through the archive and adjust ranks based on newest insertion.
******************************************************************************/
void BEERS::AdjustRanks(ArchiveStruct * pCur, ArchiveStruct * pArch)
{
   ArchiveStruct * p1;

   for(p1 = pArch; p1 != NULL; p1 = p1->pNext)
   {
      if(p1 != pCur)
      {
         if(p1->F[0] <= pCur->F[0])
         {
            (pCur->Z)++;
         }
         else
         {
            (p1->Z)++;
         }
      }
   }/* end for() */
}/* end AdjustRanks() */

/******************************************************************************
AssignModelProbs()

Compute the likelihood of each model in the current archive.
******************************************************************************/
void BEERS::AssignModelProbs(ArchiveStruct * pA, double Fbest)
{
   ArchiveStruct * p1;

   for(p1 = pA; p1 != NULL; p1 = p1->pNext)
   {
      if(Fbest <= 0.00)
      {
         p1->P = 1.00/(p1->Z + 1);
      }
      else
      {
         p1->P = Fbest/(p1->F[0]);
      }
   }/* end for() */
}/* end AssignModelProbs() */

/******************************************************************************
CalcProbabilities()

Calculate the probability of accepting a candidate solution based on 
explorative and exploitative behavior.
******************************************************************************/
void BEERS::CalcProbabilities(ArchiveStruct * pC, ArchiveStruct * pA, double dmax, double * pMin, double * pRange, double * pExploit, double * pExplore)
{
   ArchiveStruct * p1;
   int np = pC->nX;
   double * x1, * x2;
   double x1i, x2i, dtst;

   *pExplore = 1.00;
   *pExploit = 1.00;
   if(m_CurSample < 2) return;

   //compute distance to nearest point in archive
   double dmin = NEARLY_HUGE;

   x2 = pC->X;
   for(p1 = pA; p1 != NULL; p1 = p1->pNext)
   {
      x1 = p1->X;

      //compute normalized distance between points
      dtst = 0.00;
      for(int i = 0; i < np; i++)
      {
         x1i = (x1[i] - pMin[i])/pRange[i];
         x2i = (x2[i] - pMin[i])/pRange[i];

         dtst += ((x2i - x1i)*(x2i - x1i));
      }/* end for() */
      dtst = sqrt(dtst);

      //re-assign min distance, if needed
      if(dtst < dmin)
      {
         dmin = dtst;
         *pExploit = p1->P;
      }/* end if() */
   }/* end for() */   

   //exploration acceptance probability is ratio of min and max distance
   *pExplore = (dmin/dmax);
}/* end CalcProbabilities() */

/******************************************************************************
EstimateMaxDistance()

Estimate the maximum distance between a candidate point and the nearest point 
in the archive.
******************************************************************************/
double BEERS::EstimateMaxDistance(ArchiveStruct * pA, double * pMin, double * pRange)
{
   ArchiveStruct * p1, * p2;
   double * x1, * x2;   
   int np;
   double x1i, x2i, dtst;
   double dmax = 0.00;

   if(pA == NULL) return 1.00;

   //compute distances between existing points
   np = pA->nX;
   for(p1 = pA; p1 != NULL; p1 = p1->pNext)
   {
      x1 = p1->X;
      for(p2 = p1->pNext; p2 != NULL; p2 = p2->pNext)
      {         
         x2 = p2->X;

         //compute normalized distance between points
         dtst = 0.00;
         for(int i = 0; i < np; i++)
         {
            x1i = (x1[i] - pMin[i])/pRange[i];
            x2i = (x2[i] - pMin[i])/pRange[i];

            dtst += ((x2i - x1i)*(x2i - x1i));
         }/* end for() */
         dtst = sqrt(dtst);

         //re-assign max distance, if needed
         if(dtst > dmax)
         {
            dmax = dtst;
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //randomly sample 1000*np points and revise max distance if needed
   int npoints = np*1000;
   for(int pt = 0; pt < npoints; pt++)
   {      
      for(p1 = pA; p1 != NULL; p1 = p1->pNext)
      {
         x1 = p1->X;

         //compute normalized distance between points
         dtst = 0.00;
         for(int i = 0; i < np; i++)
         {
            x1i = (x1[i] - pMin[i])/pRange[i];
            x2i = ((double)MyRand() / (double)MY_RAND_MAX);

            dtst += ((x2i - x1i)*(x2i - x1i));
         }/* end for() */
         dtst = sqrt(dtst);

         //re-assign max distance, if needed
         if(dtst > dmax)
         {
            dmax = dtst;
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //printf("estimated dmax = %E\n", dmax);
   
   return dmax; 
}/* end EstimateMaxDistance() */

/******************************************************************************
WriteArchive()

Write out archive.
******************************************************************************/
void BEERS::WriteArchive(void) 
{
   FILE * pOut = fopen("OstArchive.txt", "w");
   for(ArchiveStruct * p1 = m_pArchive; p1 != NULL; p1 = p1->pNext)
   {
      fprintf(pOut, "%E  ", p1->F[0]);
      fprintf(pOut, "%04d  ", (int)(p1->Z));
      fprintf(pOut, "%E  ", p1->P);
      for(int i = 0; i < p1->nX; i++)
      {
         fprintf(pOut, "%E  ", p1->X[i]);
      }
      fprintf(pOut, "\n");
   }
   fclose(pOut);
}/* end WriteArchive() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void BEERS::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Balanced Exploration-Exploitation Random Search\n");
   fprintf(pFile, "Desired Samples         : %d\n", m_NumSamples);
   fprintf(pFile, "Actual Samples          : %d\n", m_CurSample);
   
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void BEERS::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];

   m_NumSamples = 25;

   //read in BEERS configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open BEERS config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginBEERS", pFileName) == true)
   {
      FindToken(pFile, "EndBEERS", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginBEERS", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndBEERS") == NULL)
      {         
         if(strstr(line, "NumSamples") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumSamples); 
         }/*end else if() */         
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
DestroyArchive()

Free up memory of archive.
******************************************************************************/
void BEERS::DestroyArchive(ArchiveStruct * pArch)
{
   ArchiveStruct * pDel;
   for(ArchiveStruct * pCur = pArch; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }
}/* end DestroyArchive() */

/******************************************************************************
BEERS_Program()

Calibrate or optimize the model using BEERS.
******************************************************************************/
void BEERS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("BEERS", 1);
   BEERS * TheAlg = new BEERS(model);
   MEM_CHECK(TheAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end BEERS_Program() */

