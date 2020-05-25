/******************************************************************************
File     : PDDS.h
Author   : L. Shawn Matott
Copyright: 2014, L. Shawn Matott

The PDDS algorithm is a parallel version of the DDS algorithm. It has been
ported from an earlier fortran implementation.

Version History
12-23-14    lsm   created file
******************************************************************************/
#ifndef PDDS_ALGORITHM_H
#define PDDS_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;
class StatsClass;

/******************************************************************************
class PDDSAlgorithm
******************************************************************************/
class PDDSAlgorithm : public AlgorithmABC
{
   public:
      PDDSAlgorithm(ModelABC * pModel);
      ~PDDSAlgorithm(void){ DBG_PRINT("PDDSAlgorithm::DTOR"); Destroy(); }

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
      void InitDdsDataMembers(void);
      void DestroyDdsDataMembers(void);
      double neigh_value(double x_cur, double x_min, double x_max, double r);
      double random_number(void);
      double obj_func(int nopt, double * x_values);
      int max_int(int a, int b);
      void MakeParameterCorrections(double * x, double * xb, int n, double a);

      ModelABC   *m_pModel; //Pointer to model being optimized
      StatsClass *m_pStats; //Pointer to statistics

      double m_r_val; //perturbation number 0<r<1
      int m_CurIter;
      int m_MaxIter; //maximum number of iterations                                                                                            
      int m_UserSeed; //random number generator seed
      bool m_UserSuppliedInit; //if true, then algorithm starts with users best guess (param->EstVal)
									    //if false, random parameter set chosen		
      int m_NumInit;
      double ** m_pInit;

      /* ==================================================================
      Store global variables from MOD_DDS.f90 as private data members.
      ================================================================== */
      int m_rank;
      int m_nprocessors;
      int m_master;

      //debug settings
      bool m_DEBUG_dds; //read from PDDS section of input file
      bool m_DEBUG_neigh_value; //read from PDDS section of input file

      // DDS user input values
      int m_num_dec;//num_dec = number of decision variables
      int m_ngd;
      int m_ign;
      char m_use_opt[400];//read from PDDS section of input file
      double m_alpha;//read from PDDS section of input file
      double m_beta;//read from PDDS section of input file

      char ** m_DVnames; //design variable names
      double * m_s_min;
      double * m_s_max;
      double * m_harvest;

      // DDS Output declarations
      double m_Fbest;
      double m_to_max; //always +1 in Ostrich
      double * m_sbest;
      double * m_stest;    
}; /* end class DDSAlgorithm */

extern "C" {
void PDDS_Program(int argC, StringType argV[]);
}

#endif /* PDDS_ALGORITHM_H */
