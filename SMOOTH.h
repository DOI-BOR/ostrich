/******************************************************************************
File     : SMOOTH.h
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

Simple Multi-Objective Optimization Test Heuristic (SMOOTH).

A simple random search algorithm for multi-objective optimization. It is useful
for testing OSTRICH's underlying MO support structures.

Version History
05-30-14    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef SMOOTH_H
#define SMOOTH_H

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;

/******************************************************************************
class SMOOTH

******************************************************************************/
class SMOOTH : public AlgorithmABC
{
   public:
      SMOOTH(ModelABC * pModel);
      ~SMOOTH(void){ DBG_PRINT("SMOOTH::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void){ return;}
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      int UpdateArchive(double * pX, int nX, double * pF, int nF);
      void EvaluateSamples(int nSamples);

      ModelABC * m_pModel;
      ArchiveStruct * m_pNonDom; //non-dominated solutions
      ArchiveStruct * m_pDom; //dominated solutions
      int m_NumNonDom; //number of non-dominated solutions
      int m_NumDom; //number of dominated solutions
      int m_SamplesPerIter;
      int m_MaxIters;
      int m_CurIter;
}; /* end class SMOOTH */

extern "C" {
void SMOOTH_Program(int argC, StringType argV[]);
}

#endif /* SMOOTH_H */


