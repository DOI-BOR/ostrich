/******************************************************************************
File      : LatinHypercube.h
Author    : L. Shawn Matott
Copyright : 2006, L. Shawn Matott

Encapsulates a lating hypercube sampling strategy for initializing populations.

Version History
06-14-06    lsm   created
******************************************************************************/
#ifndef LATIN_HYPERCUBE_H
#define LATIN_HYPERCUBE_H

#include "MyHeaderInc.h"

/******************************************************************************
class LatinHypercube

Encapsulates a Latin Hypercube Sampling strategy.
******************************************************************************/
class LatinHypercube
{
   public:      
      LatinHypercube(int rows, int cols);
      ~LatinHypercube(void){ DBG_PRINT("LatinHypercube::DTOR"); Destroy(); }
      void Destroy(void);   
      void InitRow(int row, double min, double max); //uniform sampling
      void InitRow(int row, double min, double max, double sd); //Gaussian sampling
      double SampleRow(int row);
      void ReDim(int cols);

   private:
      double ** m_pVals;
      int    *  m_pCount;
      int m_Rows;
      int m_Cols;
      int m_MaxCols;
}; /* end class LatinHypercube */

#endif /* LATIN_HYPERCUBE_H */

