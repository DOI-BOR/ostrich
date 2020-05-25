/******************************************************************************
File      : ComboSA.h
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

An implementation of the simulated annealing algorithm for use with combinatorial
and integer problems.

Version History
10-19-05    lsm   created
01-01-07    lsm   Algorithm uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef COMBO_SA_H
#define COMBO_SA_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward decs
class ModelABC;
class ModelBackup;
class ParameterABC;
class ParameterGroup;
class StatsClass;

/******************************************************************************
class ComboSA

******************************************************************************/
class ComboSA : public AlgorithmABC
{
   public:
      ComboSA(ModelABC * pModel);
      ~ComboSA(void){ DBG_PRINT("ComboSA::DTOR"); Destroy(); }
      void Destroy(void);

      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_NumOuter; }

   private :
      void GenerateRandomMove(ParameterABC * pParam);
      void GenerateRandomMove(void);
      double Transition(double initVal);
      double Equilibrate(double initVal);
      double Melt(double initVal);
      void StoreBest(void);
      void RestoreBest(void);

      int m_NumOuter; //current iteration number
      int m_MaxOuter; //maximum outer iterations
      int m_MaxInner; //maximum inner iterations
      double m_InitTemp; //Initial temperature
      double m_CurTemp;  //Current temperature
      double m_TempFactor; //Temperature reduction factor
      double m_dEavg;      //average energy change, estimated from melting phase
      double m_StopVal;  //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)
      ModelABC * m_pModel;
      ModelBackup * m_pTransBackup;
      double * m_pMelts;
      double * m_Finner;

      int m_NumMelts; //num. of obj. func. evals. used to determine initial temperature

      double * m_pBest;/* array containg current best parameter set */
      
      StatsClass * m_pStats;

      //metrics
      int m_MeltCount;
      int m_TransCount;
      int m_NumAborts;
      int m_EquilCount;
      int m_NumUprViols;
      int m_NumLwrViols;
      int m_NumUphill;
      int m_NumDownhill;
      double m_CurProb;
      double m_InitProb;
      double m_TotProb;
      int m_NumProbTests;
}; /* end class ComboSA */

extern "C" {
void CSA_Program(int argc, StringType argv[]);
}

#endif /* SA_ALGORITHM_H */

