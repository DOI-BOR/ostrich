/******************************************************************************
File     : GLUE.h
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Generalized Likelihood Uncertainty Engine - GLUE

Version History
06-23-10    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef GLUE_H
#define GLUE_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class ModelABC;
class ParameterABC;

/******************************************************************************
class GLUE

******************************************************************************/
class GLUE : public AlgorithmABC
{
   public:
      GLUE(ModelABC * pModel);
      ~GLUE(void){ DBG_PRINT("GLUE::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      
      void WarmStart(void){ return;}
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      void EvaluateSamples(void);
      //void BcastSamples(void);
      void EvalSamplesParallel(void);

      ModelABC * m_pModel;
      SampleStruct * m_pSamples;
      SampleStruct * m_pBehavioral;
      long long m_MaxSamples;
      int m_NumDesired;
      int m_NumFound;
      int m_SamplesPerIter;
      int m_CurIter;
      double m_Threshold;
      //loop counters for sample processing
      int m_iStart; 
      int m_iEnd;
      int * m_iCounts;
      int * m_iDispls;

      //buffers used in MPI-parallel communication
      //double * m_pBuf;
      double * m_pMyBuf;
      //double * m_pTmpBuf;
      double * m_pBigBuf;
}; /* end class GLUE */

extern "C" {
void GLUE_Program(int argC, StringType argV[]);
}

#endif /* GLUE_H */


