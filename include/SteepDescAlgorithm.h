/******************************************************************************
File      : SteepDescAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of the Steepest Descent optimization algorithm.

Version History
08-26-03    lsm   created 
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef STEEP_DESC_ALGORITHM_H
#define STEEP_DESC_ALGORITHM_H

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;
class StatsClass;
class OptSearchClass;
class OptMathClass;

/******************************************************************************
class SteepDescAlgorithm

The Steepest Descent algorithm is a first-order optimization algorithm, which 
utilizes the negative of the gradient to determine the optimial search 
direction. 
******************************************************************************/
class SteepDescAlgorithm : public AlgorithmABC
{
   public:
      SteepDescAlgorithm(ModelABC * pModel);
      ~SteepDescAlgorithm(void){ DBG_PRINT("SteepDescAlgorithm::DTOR"); Destroy(); }
      void Destroy(void);
      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      //max. # of iterations, where one 1D search is counted as an iteration
      int m_MaxIter; 

      /* if difference between obj. function from the previous iteration 
      is less than the convergence value, the algorithm exits.*/
      double m_ConvVal; 

      int m_NumParams; //numer of parameters/design vars

      double * m_pSearchDir; //search direction

      ModelABC * m_pModel;
      StatsClass * m_pStats; //calibration statistics
      OptMathClass * m_pMath;
      OptSearchClass * m_pSearchAlg; //1-dimensional search algorithm

      //metrics
      int m_AlgCount;
      int m_CurIter;
      int m_NumUprViols;
      int m_NumLwrViols;
}; /* end class PowellAlgorithm */

extern "C" {
void STPDSC_Program(int argc, StringType argv[]);
}

#endif /* STEEP_DESC_ALGORITHM_H */

