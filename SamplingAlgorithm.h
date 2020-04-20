/******************************************************************************
File      : SamplingAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2007, L. Shawn Matott

An implementation of a sampling algorithm. Loosely based on Big Bang-Big Crunch
(BB-BC).

Version History
11-21-07    lsm   created 
******************************************************************************/
#ifndef SAMPLING_ALGORITHM_H
#define SAMPLING_ALGORITHM_H

#include "MyHeaderInc.h"

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;
class StatsClass;

/******************************************************************************
class SamplingAlgorithm

An implementation of a sampling algorithm. Loosely based on Big Bang-Big Crunch
(BB-BC).
******************************************************************************/
class SamplingAlgorithm : public AlgorithmABC
{
   public:
      SamplingAlgorithm(ModelABC * pModel);
      ~SamplingAlgorithm(void){ DBG_PRINT("SamplingAlgorithm::DTOR"); Destroy(); }
      void Destroy(void);
      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      void UserDefinedEvaluations(void);
      void GenerateInitialSamples(int num);
      ParameterList * InsertParamSet(double F);      

      /* ---------------------------------------
      if true, initialize using LHS sample. 
      else, use best user-supplied configuration
      --------------------------------------- */
      bool m_bRndInit;  

      double m_Radius;  //search radius
      double *  m_pSD;  //parameter std. deviations
      double *  m_pFwd; //parameter forward perturbation weights (prob. that dX is positive)
      double *  m_pLwr; //Xmin
      double *  m_pUpr; //Xmax

      int m_MaxEvals;   //targeted maximum number of model evaluations
      int m_MaxIter;    //maximum number of iterations [ sqrt(m_MaxEvals)]
      int m_NumSamples; //number samples per iteration [ m_MaxEvals/m_MaxIter ]
      int m_NumExtra;   //number of extra initial samples [ m_MaxEvals - (m_MaxIter*m_NumSamples) ]   

      int m_NumParams;

      ParameterList * m_pAll; //list of all parameter sets
      MyPoint * m_pSamples;   //list of new samples to be evaluated
      MyPoint * m_pBest;      //pointer to best overall parameter configuration

      ModelABC * m_pModel;
      StatsClass * m_pStats; //calibration statistics

      //metrics
      int m_AlgCount;
      int m_CurIter;
}; /* end class SamplingAlgorithm */

extern "C" {
void SMP_Program(int argc, StringType argv[]);
}

#endif /* SAMPLING_ALGORITHM_H */

