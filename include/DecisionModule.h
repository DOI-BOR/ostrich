/******************************************************************************
File     : DecisionModeule.h
Author   : L. Shawn Matott
Copyright: 2006, L. Shawn Matott

When a surrogate-based approach is used, the Decision Module determines which
of the set of models of varying complexity should be executed.

Version History
04-04-06    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef DECISION_MODULE_H
#define DECISION_MODULE_H

#include "MyHeaderInc.h"

//forward declarations
class ModelABC;

/******************************************************************************
class DecisionModule

Stores all possible models and a database of all model evaluations.
*****************************************************************************/
class DecisionModule
{   
   public:
      DecisionModule(ModelABC * pComplex);
      ~DecisionModule(void){ DBG_PRINT("DecisionModule::DTOR"); Destroy(); }
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      double Execute(void);
      void Bookkeep(bool bFinal);

   private:
      int GetRank(int i, double * A, int size, int type);
      int GetRank(int i, int * A, int size, int type);
      void CollectMetrics(void);
      void ShareWeights(void);
      double EvalComplex(void);
      double EvalSurrogate(int model_id);
      int GetNumModelParams(int id);

      ModelABC *  m_pComplex; //complex model
      ModelABC ** m_pModels;  //surrogate models
      int * m_pWeights; //selection weights

      int m_NumModels; //number of models (complex + surrogates)

      int m_SelectionScheme; //type of selection scheme, either biased or unbiased

      //metrics
      int m_TotalEvals;
      int * m_pEvals; //number of evals for each model
      double * m_pBestAICc; //best AICc for each model
      double * m_pBestWSSE; //best WSSE for each model
      double * m_pAICc; //temporary storage for AICc for each model
      double * m_pWSSE; //temporary storage for WSSE for each model
}; /* end class DecisionModule */

#endif /* DECISION_MODULE_H */


