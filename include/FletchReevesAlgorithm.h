/******************************************************************************
File      : FletchReevesAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of the Fletcher-Reeves optimization algorithm.

Version History
08-26-03    lsm   created 
02-17-04    lsm   switched implementation
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC)
******************************************************************************/
#ifndef FLETCH_REEVES_ALGORITHM_H
#define FLETCH_REEVES_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class StatsClass;
class ModelABC;
class OptMathClass;
class OptSearchClass;

/******************************************************************************
class FletchReevesAlgorithm

The Fletcher-Reeves algorithm is a first-order optimization algorithm, which 
utilizes the concept of conjugate directions in conjunction with the steepest-
descent information (negative of the gradient). 
******************************************************************************/
class FletchReevesAlgorithm : public AlgorithmABC
{
   public:
      FletchReevesAlgorithm(ModelABC * pModel);
      ~FletchReevesAlgorithm(void){ DBG_PRINT("FletchReevesAlgorithm::DTOR"); Destroy(); }
      void Destroy(void);

      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);

      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      //max. # of iterations, where one 1D search is counted as an iteration
      int m_MaxIter;
      int m_MaxCount;

      int m_CurIter; 

      /* if difference between obj. function from the previous iteration 
      is less than the convergence value, the algorithm exits.*/
      double m_ConvVal; 

      int m_NumParams; //numer of parameters/design vars

      ModelABC * m_pModel;
      StatsClass * m_pStats; //calibration statistics
      OptMathClass * m_pMath;
      OptSearchClass * m_pSearchAlg; //1-dimensional search algorithm

      //metrics
      int m_NumRestarts;
      int m_NumUprViols;
      int m_NumLwrViols;
      int m_AlgCount;
}; /* end class FletchReevesAlgorithm */

extern "C" {
void FLRV_Program(int argc, StringType argv[]);
}

#endif /* FLETCH_REEVES_ALGORITHM_H */

