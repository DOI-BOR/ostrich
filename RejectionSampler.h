/******************************************************************************
File     : RejectionSampler.h
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Rejection Sampling Algorithm

Version History
06-23-10    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef REJECTION_SAMPLER_H
#define REJECTION_SAMPLER_H

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;

/******************************************************************************
class RejectionSampler

******************************************************************************/
class RejectionSampler : public AlgorithmABC
{
   public:
      RejectionSampler(ModelABC * pModel, bool bMCMC);
      ~RejectionSampler(void){ DBG_PRINT("RejectionSampler::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      
      void WarmStart(void){ return;}
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      double ComputeLikelihoodRatio(double wsse);
      void EvaluateSamples(void);
      void EvalSamplesSuperMUSE(void);
      void BcastSamples(void);
      void EvalSamplesParallel(void);
      void EvalSamplesParallelAsync(void);

      ModelABC * m_pModel;
      SampleStruct * m_pSamples;
      SampleStruct * m_pAccepted;
      long long m_MaxSamples;
      int m_NumDesired;
      int m_NumBurnIn;
      int m_NumFound;
      int m_SamplesPerIter;
      int m_CurIter;
      double m_MinWSSE; //for rejection sampler
      double m_LastWSSE; //for Metropolis sampler
      bool m_bStedinger; //if true, use Stedinger derivation
      bool m_bBeven; //if true, use Beven derivation
      double m_ShapeFactor; //shaping factor for Beven's Pseudo-Likelihood
      bool m_bMetropolis; //if true, use Metropolis MCMC
      double m_TelescopeRate;

      //buffers used in MPI-parallel communication
      double * m_pBuf;
      double * m_pMyBuf;
      double * m_pTmpBuf;
      double * m_pBigBuf;
}; /* end class RejectionSampler */

extern "C" {
void RJSMP_Program(int argC, StringType argV[]);
void METRO_Program(int argC, StringType argV[]);
}

#endif /* REJECTION_SAMPLER_H */


