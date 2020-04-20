/******************************************************************************
File      : PowellAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of Powell's optimization algorithm.

Version History
08-26-03    lsm   created 
02-17-04    lsm   switched implementation
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC

******************************************************************************/
#ifndef POWELL_ALGORITHM_H
#define POWELL_ALGORITHM_H

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;
class StatsClass;
class OptSearchClass;

/******************************************************************************
class PowellAlgorithm

Powell's method is a zero-order optimization algorithm, which utilizes the
concept of conjugate directions to determine the optimial search direction.
******************************************************************************/
class PowellAlgorithm : public AlgorithmABC
{
   public:
      PowellAlgorithm(ModelABC * pModel);
      ~PowellAlgorithm(void){ DBG_PRINT("PowellAlgorithm::DTOR"); Destroy(); }
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

      int m_NumDirs; //numer of search directions

      double ** m_pSearchDirs; //array of search directions

      ModelABC * m_pModel;
      StatsClass * m_pStats; //calibration statistics
      OptSearchClass * m_pSearch; //1-dimensional search

      //metrics
      int m_AlgCount;
      int m_NumRestarts;
      int m_NumUprViols;
      int m_NumLwrViols;
      int m_CurIter;

}; /* end class PowellAlgorithm */

extern "C" {
void PWL_Program(int argc, StringType argv[]);
}

#endif /* POWELL_ALGORITHM_H */

