/******************************************************************************
File     : DDS.h
Author   : James R. Craig and Bryan Tolson
Copyright: 2006, James R. Craig and Bryan Tolson

The DDS algorithm uses an intelligent greedy search to search the parameter 
space and find the global minimum of an optimization problem. The algorithm is 
discussed in Tolson and Shoemaker, Water Resources Research, 2007

Dynamically dimensioned Search (DDS) version 1.1 algorithm by Bryan Tolson
c++ version (original was coded in Matlab, translated to Fortran by Bryan Tolson)
Translated to c++ by James Craig (July-Aug 2006)
DDS is an n-dimensional continuous global optimization algorithm

Ostrich updates required in following files:
	MyTypes.h
	Ostrich.cpp
	Utilities.h
	Utilities.cpp

Version History
03-01-06    jrc   created file
******************************************************************************/
#ifndef DDS_ALGORITHM_H
#define DDS_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class StatsClass;
class ModelABC;
class ParameterABC;

/******************************************************************************
class DDSAlgorithm
******************************************************************************/
class DDSAlgorithm : public AlgorithmABC
{
   public:
      DDSAlgorithm(ModelABC * pModel);
      ~DDSAlgorithm(void){ DBG_PRINT("DDSAlgorithm::DTOR"); Destroy(); }

      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);

      void SetPerturbationValue(double val) { m_r_val = val; }
      void ResetUserSeed(int seed);
      void SetNoUserInit(void){ m_UserSuppliedInit = false; }
      void SetBudget(int budget){ m_MaxIter = budget; }

      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      ModelABC   *m_pModel;             //Pointer to model being optimized
		StatsClass *m_pStats;             //Pointer to statistics

		double m_r_val;							//perturbation number 0<r<1
		int m_MaxIter;						//maximum number of iterations                                                                                            
      int m_CurIter;
		int m_UserSeed;					  //random number generator seed
		bool m_UserSuppliedInit;		//if true, then algorithm starts with users best guess (param->EstVal)
																				//if false, random parameter set chosen		

			double PerturbParam(const double &best_value, ParameterABC * pParam);  
         void MakeParameterCorrections(double * x, double * xb, int n, double a);

}; /* end class DDSAlgorithm */

extern "C" {
void DDS_Program(int argC, StringType argV[]);
}

#endif /* DDS_ALGORITHM_H */
