/******************************************************************************
File     : DDSAU.h
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

This algorithm seeks to identify behavioral parameters set by repeatedly applying 
a DDS search from alternative starting points in the parameter space. An optional 
group will configure the DDS for Approximation of Uncertainty (DDS-AU) algorithm 
and will be processed if ProgramType is set to “DDSAU”.

Version History
11-15-15    lsm   created file
******************************************************************************/
#ifndef DDSAU_ALGORITHM_H
#define DDSAU_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class StatsClass;
class ModelABC;
class ParameterABC;

/******************************************************************************
class DDSAU
******************************************************************************/
class DDSAU : public AlgorithmABC
{
   public:
      DDSAU(ModelABC * pModel);
      ~DDSAU(void){ DBG_PRINT("DDSAU::DTOR"); Destroy(); }

      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void){ return;}

      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      void OptimizeSerial(void);
      void OptimizeParallel(void);

      ModelABC   *m_pModel;             //Pointer to model being optimized
		StatsClass *m_pStats;             //Pointer to statistics
     
      char ** m_pBehavioral;
      double * m_fBehavioral;
 
		double m_r_val; //perturbation number 0<r<1
      int m_nsols; //number of searches
      int m_nbhvr; //number of behavioral smaples
      int m_MinIter; //minimum number of iterations per search
		int m_MaxIter; //maximum number of iterations per search
      int m_CurIter;
      bool m_bParallel; //true == perform searches in parallel
      double m_fmax; //fitness threshold
      bool m_bRandomize; //true == randomize the behavioral solutions
      bool m_bReviseAU; // true -- revise previous search

		double PerturbParam(double * best_value, ParameterABC * pParam);  
      void MakeParameterCorrections(double * x, double * xb, int n, double a);

}; /* end class DDSAU */

extern "C" {
void DDSAU_Program(int argC, StringType argV[]);
}

#endif /* DDSAU_ALGORITHM_H */
