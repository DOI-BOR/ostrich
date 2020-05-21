/******************************************************************************
File      : BisectionAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2007, L. Shawn Matott

An implementation of a basic bisection algorithm.

Version History
08-01-07    lsm   created 
******************************************************************************/
#ifndef BISECTION_ALGORITHM_H
#define BISECTION_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class ModelABC;
class ParameterGroup;
class StatsClass;

/******************************************************************************
class BisectionAlgorithm

The bisection algorithm optimizes one parameter at a time via simple bisection.
******************************************************************************/
class BisectionAlgorithm : public AlgorithmABC
{
   public:
      BisectionAlgorithm(ModelABC * pModel);
      ~BisectionAlgorithm(void){ DBG_PRINT("BisectionAlgorithm::DTOR"); Destroy(); }
      void Destroy(void);
      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      int m_MaxOuter; 
      int m_MaxInner; 

      int m_NumParams;

      /* if difference between obj. function from the previous iteration 
      is less than the convergence value, the algorithm exits.*/
      double m_ConvVal; 

      ModelABC * m_pModel;
      StatsClass * m_pStats; //calibration statistics

      //metrics
      int m_AlgCount;
      int m_CurIter;
}; /* end class PowellAlgorithm */

extern "C" {
void BIS_Program(int argc, StringType argv[]);
}

#endif /* BISECTION_ALGORITHM_H */

