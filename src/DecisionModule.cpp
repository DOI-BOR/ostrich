/******************************************************************************
File     : DecisionModeule.h
Author   : L. Shawn Matott
Copyright: 2006, L. Shawn Matott

When a surrogate-based approach is used, the Decision Module determines which
of the set of models of varying complexity should be executed.

Version History
04-04-06    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "DecisionModule.h"
#include "Model.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "SurrogateParameterGroup.h"
#include "ParameterGroup.h"

#include "Utility.h"
#include "Exception.h"

/* ---------------------------------------------------
Two different model selection schemes:
   UNBIASED : no randomness, evaluate each model before 
   making selection, correct decision is guranteed but
   method is computationally epxensive.

   BIASED_RANDOM : Biased-but-random selection, only
   one model is selected using an adaptive biased-but-
   random weighting scheme
 --------------------------------------------------- */
#define SELECTION_SCHEME_UNBIASED        (0)
#define SELECTION_SCHEME_BIASED_RANDOM   (1)

#define RANK_TYPE_DESCENDING   (0)
#define RANK_TYPE_ASCENDING    (1)

/******************************************************************************
CTOR

Constructs a DecisionModule by reading the SurrogateModels section of the 
input file.
******************************************************************************/
DecisionModule::DecisionModule(ModelABC * pComplex)
{
   int i;
   FILE * pFile;
   char * line, tmpStr[DEF_STR_SZ], modeStr[DEF_STR_SZ];
   char beg_tok[DEF_STR_SZ], end_tok[DEF_STR_SZ];

   m_SelectionScheme = SELECTION_SCHEME_BIASED_RANDOM;
   //m_pWeights = NULL;

   m_pComplex = pComplex;
   m_pModels = NULL;
   m_NumModels = 1;
   m_TotalEvals = 0;
   m_pEvals = NULL;
   m_pBestAICc = NULL;
   m_pBestWSSE = NULL;
   m_pAICc = NULL;
   m_pWSSE = NULL;
  
   pFile = fopen(GetSrgFileName(), "r");
   
   if(CheckToken(pFile, "SelectionScheme", GetSrgFileName()) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmpStr, modeStr);
      MyStrLwr(modeStr);
      if(strcmp(modeStr, "unbiased") == 0){
         m_SelectionScheme = SELECTION_SCHEME_UNBIASED;
         printf("Selection scheme is unbiased, all models evaluated every time\n");}
      else if(strcmp(modeStr, "biased-but-random") == 0){
         m_SelectionScheme = SELECTION_SCHEME_BIASED_RANDOM;
	      printf("Selection scheme is biased-but-random\n");}
      else {
         m_SelectionScheme = SELECTION_SCHEME_BIASED_RANDOM;
	      printf("Unknown selection scheme: %s\n", modeStr);
         printf("defaulting to biased-but-random.\n");}
   }
   rewind(pFile);
   if(CheckToken(pFile, "NumberOfSurrogates", GetSrgFileName()) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %d", tmpStr, &m_NumModels);
      m_NumModels++;
   }
   rewind(pFile);

   for(i = 1; i < m_NumModels; i++)
   {
      /* -------------------------------------------------
      assemble dynamic tokens and use FindToken() to make 
      sure all surrogate sections are in the input file
      ------------------------------------------------- */
      sprintf(beg_tok, "Begin_S%d_Model", i);
      sprintf(end_tok, "End_S%d_Model", i);
      FindToken(pFile, beg_tok, GetSrgFileName());
      FindToken(pFile, end_tok, GetSrgFileName());
      rewind(pFile);
   }
   fclose(pFile);

   m_pAICc = new double[m_NumModels]; //temporary storage of AICc values
   m_pWSSE = new double[m_NumModels]; //temporary storage of WSSE values
   m_pWeights = new int[m_NumModels]; //array of selection weights

   //create the dynamically-named surrogate models
   m_pModels = new ModelABC *[m_NumModels];
   MEM_CHECK(m_pModels);

   m_pModels[0] = m_pComplex;
   for(i = 1; i < m_NumModels; i++)
   {
      sprintf(tmpStr, "S%d", i);
      m_pModels[i] = new SurrogateModel(GetDynFileName(tmpStr), m_pComplex, tmpStr);
      MEM_CHECK(m_pModels[i]);
   }
   GetDynFileName(NULL);//clean up last dynamic file

   m_pEvals = new int[m_NumModels];
   m_pBestWSSE = new double[m_NumModels];
   m_pBestAICc = new double[m_NumModels];
   MEM_CHECK(m_pBestAICc);

   for(i = 0; i < m_NumModels; i++)
   {
      m_pBestAICc[i] = NEARLY_HUGE;
      m_pBestWSSE[i] = NEARLY_HUGE;
      m_pEvals[i] = 0;
   }

   //delete temporary input file
   remove(GetSrgFileName());

   IncCtorCount();
} /* end DecisionModule::CTOR */

/******************************************************************************
DTOR

Destroys a DecisionModule by freeing up memory.
******************************************************************************/
void DecisionModule::Destroy(void)
{
   int i;
   for(i = 1; i < m_NumModels; i++)
   {
      delete m_pModels[i];
   }
   delete [] m_pEvals;
   delete [] m_pBestAICc;
   delete [] m_pBestWSSE;
   delete [] m_pAICc;
   delete [] m_pWSSE;
   delete [] m_pWeights;
   IncDtorCount();
} /* end DecisionModule::DTOR */

/******************************************************************************
WriteMetrics()

Write out metrics.
******************************************************************************/
void DecisionModule::WriteMetrics(FILE * pFile)
{
   int i;

   if(m_SelectionScheme == SELECTION_SCHEME_UNBIASED)
   {
      fprintf(pFile, "Selection Scheme       : unbiased\n");
   }
   else
   {
      fprintf(pFile, "Selection Scheme       : biased-but-random\n");
      fprintf(pFile, "Complex Selection Weight   : %d\n", m_pWeights[0]);
      for(i = 1; i < m_NumModels; i++)
      {
         fprintf(pFile, "S%02d Selection Weight   : %d\n",i, m_pWeights[i]);
      }
   }
   fprintf(pFile, "Complex Evals           : %d\n", m_pEvals[0]);
   fprintf(pFile, "Complex Best AICc       : %lf\n", m_pBestAICc[0]);
   fprintf(pFile, "Complex Best WSSE       : %lf\n", m_pBestWSSE[0]);

   for(i = 1; i < m_NumModels; i++)
   {
      fprintf(pFile, "S%02d Evals               : %d\n", i, m_pEvals[i]);
      fprintf(pFile, "S%02d Best AICc           : %lf\n", i, m_pBestAICc[i]);
      fprintf(pFile, "S%02d Best WSSE           : %lf\n", i, m_pBestWSSE[i]);
   }/* end for() */

   fprintf(pFile, "Total Evals             : %d\n", m_TotalEvals);
}/* end WriteMetrics() */

/******************************************************************************
GetRank()

Compute the rank of element i in array A.
******************************************************************************/
int DecisionModule::GetRank(int i, int * A, int size, int type)
{
   int j;
   int rank = 1;

   for(j = 0; j < size; j++)
   {
      if((A[j] < A[i]) && (type == RANK_TYPE_ASCENDING)) { rank++;}
      if((A[j] > A[i]) && (type == RANK_TYPE_DESCENDING)){ rank++;}
   }/* end for() */
   return rank;
}/* end GetRank() */

/******************************************************************************
GetRank()

Compute the rank of element i in array A.
******************************************************************************/
int DecisionModule::GetRank(int i, double * A, int size, int type)
{
   int j;
   int rank = 1;

   for(j = 0; j < size; j++)
   {
      if((A[j] < A[i]) && (type == RANK_TYPE_ASCENDING)) { rank++;}
      if((A[j] > A[i]) && (type == RANK_TYPE_DESCENDING)){ rank++;}
   }/* end for() */
   return rank;
}/* end GetRank() */

/******************************************************************************
Execute()

Run the appropriate model.
******************************************************************************/
double DecisionModule::Execute(void)
{
   static bool bFirstTime = true;
   double AICc, WSSE, minAICc;
   double F;   
   int nobs, npi, id, i, sum;
   int model_id, rand;

   nobs = m_pComplex->GetObsGroupPtr()->GetNumObs();

   /* -------------------------------------------------------------------------
   Assign initial selection weights or apply unbiased selection.
   ------------------------------------------------------------------------- */
   if((bFirstTime == true) || (m_SelectionScheme == SELECTION_SCHEME_UNBIASED))
   {
      /* --------------------------------------------------
      Evaluate each model
      --------------------------------------------------- */
      for(model_id = 0; model_id < m_NumModels; model_id++)
      {
         npi = GetNumModelParams(model_id) + 1;
         if(model_id > 0) WSSE = m_pModels[model_id]->Execute();
         else             WSSE = ((Model *)m_pComplex)->StdExecute(0.00);
         m_pEvals[model_id]++;
         m_TotalEvals++;
         AICc = ((nobs * log(WSSE/nobs)) + 2.00 * npi + 
                ((2.00 * npi * (npi + 1.00))/(nobs - npi - 1.00)));
         m_pWSSE[model_id] = WSSE;
         m_pAICc[model_id] = AICc;

         if(bFirstTime == true) //initialize 'best' WSSE and AICc
         {
            m_pBestWSSE[model_id] = WSSE;
            m_pBestAICc[model_id] = AICc;
         }
      }/* end for() */

      /* --------------------------
      Determine the best result.
      --------------------------- */
      id = 0;
      minAICc = m_pAICc[0];
      for(i = 1; i < m_NumModels; i++)
      {
         if(m_pAICc[i] <= minAICc)
         {
            id = i;
         }/* end if() */
      }/* end for() */

      bFirstTime = false;
   }/* end if() */
   else /* biased-but-random selection */
   {
      rand = MyRand() % (m_pWeights[m_NumModels]);
      sum = 0;
      id = (m_NumModels-1);

      for(i = 0; i < (m_NumModels-1); i++)
      {
         sum += m_pWeights[i];

         if(rand < sum){ 
            id = i; 
            break;}
      }
   }/* end else() */
 
   /* --------------------------------
   Perform 'official' model evaluation
   -------------------------------- */
   if(id == 0) F = EvalComplex();
   else        F = EvalSurrogate(id);

   //compute AICc
   npi = GetNumModelParams(id) + 1;
   AICc = ((nobs * log(F/nobs)) + 2.00 * npi + 
         ((2.00 * npi * (npi + 1.00))/(nobs - npi - 1.00)));

   /* -----------------------------------------------------------
   Update best WSSE and AICc metrics, if applicable
   ------------------------------------------------------------ */
   if(AICc < m_pBestAICc[id]) 
   {
      m_pBestAICc[id] = AICc;
      m_pBestWSSE[id] = F;
   }

   /* ------------------------------------------------------------------------------
   Assign/revise selection weights
   ------------------------------------------------------------------------------ */
   if(m_SelectionScheme == SELECTION_SCHEME_BIASED_RANDOM)
   {
      sum = 0;
      for(i = 0; i < m_NumModels; i++)
      {
         m_pWeights[i] = GetRank(i, m_pBestAICc, m_NumModels, RANK_TYPE_DESCENDING);
         sum += m_pWeights[i];
      }
      m_pWeights[m_NumModels] = sum;
   }/* end if() */

   //bump master counter
   m_TotalEvals++;
   return F;
}/* end Execute() */

/******************************************************************************
EvalComplex()

Run the complex model and update the database.
******************************************************************************/
double DecisionModule::EvalComplex(void)
{
   double Fcmx;

   //evaluate model
   Fcmx = ((Model *)m_pComplex)->StdExecute(0.00);
   m_pEvals[0]++;

   return Fcmx;
}/* end EvalComplex() */

/******************************************************************************
EvalSurrogate()

Run the surrogate model and update the database.
******************************************************************************/
double DecisionModule::EvalSurrogate(int model_id)
{
   int i;
   int nobs;
   ModelABC * pModel;
   double Fsrg, obs;

   pModel = m_pModels[model_id];

   //execute simple model
   Fsrg = pModel->Execute();
   m_pEvals[model_id]++;

   //propagate results to complex model
   m_pComplex->SetObjFuncVal(Fsrg);

   //initialize number of observations
   nobs = pModel->GetObsGroupPtr()->GetNumObs();

   //propagate observations to complex model
   for(i = 0; i < nobs; i++)
   {
      obs = pModel->GetObsGroupPtr()->GetObsPtr(i)->GetComputedVal(false, false);
      m_pComplex->GetObsGroupPtr()->GetObsPtr(i)->SetComputedVal(obs);
   }/* end for() */

   return Fsrg;
}/* end EvalSurrogate() */

/******************************************************************************
GetNumModelParams()

Get the number of parameters in the model.
******************************************************************************/
int DecisionModule::GetNumModelParams(int id)
{
   SurrogateModel * pSrg;

   if(id == 0) return m_pComplex->GetParamGroupPtr()->GetNumParams();

   pSrg = (SurrogateModel *)(m_pModels[id]);
   return pSrg->GetSurrogateParamGroupPtr()->GetNumTiedParams();
}/* end GetNumModelParams() */

/*****************************************************************************
Bookkeep()
   Performs bookkeeping operations related to parallel executing.
******************************************************************************/
void DecisionModule::Bookkeep(bool bFinal)
{
   ShareWeights();
   if(bFinal == true)
   {
      CollectMetrics();
   }
}/* end Bookkeep() */

/*****************************************************************************
CollectMetrics()
   Collects metrics from other processors.
******************************************************************************/
void DecisionModule::CollectMetrics(void)
{     
   int id, nprocs, temp, i, j;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (nprocs == 1) return;

   /* ---------------------------------------------------------------
   Collect individual and total evals 
   --------------------------------------------------------------- */
   for(i = 1; i < nprocs; i++)
   {
      temp = m_TotalEvals;
      MPI_Bcast(&temp, 1, MPI_INTEGER, i, MPI_COMM_WORLD);
      if(id == 0) m_TotalEvals += temp;

      for(j = 0; j < m_NumModels; j++)
      {
         temp = m_pEvals[j];
         MPI_Bcast(&temp, 1, MPI_INTEGER, i, MPI_COMM_WORLD);
         if(id == 0) m_pEvals[j] += temp;
      }
   }
} /* end CollectMetrics() */

/*****************************************************************************
ShareWeights()
   When in parallel, have processors share their weights 
   (i.e. best AICc values).
******************************************************************************/
void DecisionModule::ShareWeights(void)
{  
   double temp;
   int i, j, nprocs, sum;

   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (nprocs == 1) return;
   if (m_SelectionScheme != SELECTION_SCHEME_BIASED_RANDOM) return;

   /* ---------------------------------------------------------------
   Have each processor broadcast it's best AICc and WSSE value for 
   each model.
   --------------------------------------------------------------- */
   for(i = 0; i < nprocs; i++)
   {
      for(j = 0; j < m_NumModels; j++)
      {
         temp = m_pBestAICc[j];
         MPI_Bcast(&temp, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
         if(temp < m_pBestAICc[j]){ m_pBestAICc[j] = temp;}

         temp = m_pBestWSSE[j];
         MPI_Bcast(&temp, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
         if(temp < m_pBestWSSE[j]){ m_pBestWSSE[j] = temp;}
      }/* end for() */
   }/* end for() */
   
   //assign weights
   sum = 0;
   for(i = 0; i < m_NumModels; i++)
   {
      m_pWeights[i] = GetRank(i, m_pBestAICc, m_NumModels, RANK_TYPE_DESCENDING);
      sum += m_pWeights[i];
   }
   m_pWeights[m_NumModels] = sum;
} /* end ShareWeights() */

