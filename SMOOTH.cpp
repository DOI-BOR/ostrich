/******************************************************************************
File     : SMOOTH.cpp
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

Simple Multi-Objective Optimization Test Heuristic (SMOOTH).

A simple random search algorithm for multi-objective optimization. It is useful
for testing OSTRICH's underlying MO support structures.

Version History
05-30-14    lsm   added copyright information and initial comments.
******************************************************************************/
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "SMOOTH.h"
#include "Model.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
SMOOTH::SMOOTH(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pNonDom = NULL;
   m_pDom = NULL;
   m_NumNonDom = 0;
   m_NumDom = 0;
   m_SamplesPerIter = 0;
   m_MaxIters = 0;
   m_CurIter = 0;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by SMOOTH and it's member variables.
******************************************************************************/
void SMOOTH::Destroy(void)
{
   ArchiveStruct * pCur, * pDel;

   for(pCur = m_pNonDom; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   for(pCur = m_pDom; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using SMOOTH.
******************************************************************************/
void SMOOTH::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Search for the pareto front using SMOOTH.
******************************************************************************/
void SMOOTH::Optimize(void)
{
   int num, nobj, result;
   double * X, * F, lwr, upr, range, r;
   StatusStruct pStatus;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   InitFromFile(GetInFileName());

   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();
   nobj = m_pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);

   WriteSetup(m_pModel, "SMOOTH - Simple Multi-Objective Optimization Test Heuristic");
   //write banner
   WriteBanner(m_pModel, "gen   ", "Convergence Value");

   //main optimization loop   
   for(int g = 0; g < m_MaxIters; g++)
   {
      pStatus.curIter = m_CurIter = g+1;
      if(IsQuit() == true){ break;}

      //random parameter sampling 
      WriteInnerEval(WRITE_SMP, m_SamplesPerIter, '.');
      for(int iS = 0; iS < m_SamplesPerIter; iS++)
      {
         X = new double[num];
         F = new double[nobj];
         for(int j = 0; j < num; j++) //for each parameter
         {
            pParam = pGroup->GetParamPtr(j);
            //generate a random between lower and upper bound
            lwr = pParam->GetLwrBnd();
            upr = pParam->GetUprBnd();
            range = upr - lwr;
            r = (double)MyRand() / (double)MY_RAND_MAX;
            X[j] = (r * range) + lwr;            
         }/* end for() */
         pGroup->WriteParams(X);
         m_pModel->Execute(F, nobj);
         result = UpdateArchive(X, num, F, nobj); 
         if(result == ARCHIVE_NON_DOM)
         {
            WriteInnerEval(iS+1, m_SamplesPerIter, '+');
         }/* end if() */
         else
         {
            WriteInnerEval(iS+1, m_SamplesPerIter, '-');
         }/* end else() */
      }/* end for() */
      WriteInnerEval(WRITE_ENDED, 0, '.');

      pStatus.pct = ((float)100.00*(float)(g+1))/(float)m_MaxIters;
      pStatus.numRuns = (g+1)*m_SamplesPerIter;
      WriteMultiObjRecord(m_pModel, (g+1), m_pNonDom, pStatus.pct);
      WriteStatus(&pStatus);
   }/* end for() */

   WriteMultiObjOptimal(m_pModel, m_pNonDom, m_pDom);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   //write algorithm metrics
   WriteAlgMetrics(this);
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void SMOOTH::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : SMOOTH - Simple Multi-Objective Optimization Test Heuristic\n");
   fprintf(pFile, "Max Iterations          : %d\n", m_MaxIters);
   fprintf(pFile, "Actual Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Samples per Iteration   : %d\n", m_SamplesPerIter);
   fprintf(pFile, "Non-Dominated Solutions : %d\n", m_NumNonDom);  
   fprintf(pFile, "Dominated Solutions     : %d\n", m_NumDom);     
   fprintf(pFile, "Sampling Method         : Uniform Random\n");

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
UpdateArchive()

Update the dominated and non-dominated archives with latest sample.
******************************************************************************/
int SMOOTH::UpdateArchive(double * pX, int nX, double * pF, int nF)
{
   int i;
   double Fcur, Ftst;
   ArchiveStruct * pArch, * pCur, * pPrev, * pNxt;
   bool bDominates, bIsDominated, bMarkForInsertion;
   pArch =  new ArchiveStruct;
   pArch->F = pF;
   pArch->X = pX;
   pArch->nX = nX;
   pArch->nF = nF;
   pArch->pNext = NULL;
   
   //first entry is always non-dominated
   if((m_NumDom == 0) && (m_NumNonDom == 0))
   {
      m_pDom = NULL;
      m_pNonDom = pArch;
      m_NumNonDom++;
      return ARCHIVE_NON_DOM;
   }

   //assume solution is non-dominated until we discover otherwise
   bMarkForInsertion = true;

   //compare against current list of non-dominated solutions
   ArchiveStruct * pDummy;
   for(pCur = m_pNonDom; pCur != NULL;)
   {
      //save next item since pCur->pNext may be changed during processing
      pDummy = pCur->pNext;

      //does new solution (Ftst) dominate the existing solution (Fcur)?
      bDominates = true;
      for(i = 0; i < pArch->nF; i++)
      {
         Fcur = pCur->F[i];
         Ftst = pArch->F[i];
         if(Fcur < Ftst)
         {
            bDominates = false;
            break;
         }/* end if() */
      }/* end for() */

      //is new solution (Ftst) dominated by an existing solution (Fcur)?
      if(bDominates == false)
      {
         bIsDominated = true;
         for(i = 0; i < pArch->nF; i++)
         {
            Fcur = pCur->F[i];
            Ftst = pArch->F[i];
            if(Ftst < Fcur)
            {
               bIsDominated = false;
               break;
            }/* end if() */
         }/* end for() */
      }/* end if() */

      /* -----------------------------------------------------------------------------
      Existing solution is dominated. Remove it from list of non-dominated solutions
      and mark new solution for insertion into the non-dominated list.
      ----------------------------------------------------------------------------- */
      if(bDominates == true)
      {
         //solution to be removed is at head of list.
         if(pCur == m_pNonDom)
         {
            m_pNonDom = pCur->pNext;
            m_NumNonDom--;
         }/* end if() */
         else //somewhere in middle or at end.
         {
            for(pPrev = m_pNonDom; pPrev->pNext != pCur;)
            {
               pPrev = pPrev->pNext;
            }
            pNxt = pCur->pNext;
            pPrev->pNext = pNxt; //removes item from list
            m_NumNonDom--;
         }/* end else() */

         //insert at head of dominated list
         pCur->pNext = NULL;
         pNxt = m_pDom;
         m_pDom = pCur;
         m_pDom->pNext = pNxt;
         m_NumDom++;         
      }/* end if() */
      /* -----------------------------------------------------------------------------
      New solution is dominated. Make note so that it is not inserted into the list
      of non-dominated solutions.
      ----------------------------------------------------------------------------- */
      else if(bIsDominated == true)
      {
         bMarkForInsertion = false;   
      }/* end if() */

      //advance to next item
      pCur = pDummy;
   }/* end for() */

   //insert new solution into list of non-dominated solutions?
   if(bMarkForInsertion == true)
   {
      //insert at head of non-dominated list
      pNxt = m_pNonDom;
      m_pNonDom = pArch;
      m_pNonDom->pNext = pNxt;
      m_NumNonDom++;               
      return ARCHIVE_NON_DOM;
   }/* end if() */
   else
   {
      //insert at head of dominated list
      pNxt = m_pDom;
      m_pDom = pArch;
      m_pDom->pNext = pNxt;
      m_NumDom++;               
      return ARCHIVE_DOM;
   }/* end else() */
}/* end UpdateArchive() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void SMOOTH::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   
   m_SamplesPerIter = 20;
   m_MaxIters = 50;

   //read in SMOOTH configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open SMOOTH config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginSMOOTH", pFileName) == true)
   {
      FindToken(pFile, "EndSMOOTH", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginSMOOTH", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndSMOOTH") == NULL)
      {         
         if(strstr(line, "SamplesPerIter") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_SamplesPerIter); 
         }/*end else if() */         
         else if(strstr(line, "NumIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxIters);
         }
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
SMOOTH_Program()

Calibrate or optimize the model using SMOOTH.
******************************************************************************/
void SMOOTH_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("SMOOTH", 1);
   SMOOTH * TheAlg = new SMOOTH(model);
   MEM_CHECK(TheAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end SMOOTH_Program() */

