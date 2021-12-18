/******************************************************************************
File     : ParaPADDS.h
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

The PADDS (Pareto Archived Dynamically Dimensioned Search) algorithm is a 
multi-objective version of the DDS algorithm. It has been ported from a C++ 
implementation of Mohammadamin Jahanpour (mjahanpo@uwaterloo.ca).

Version History
05-08-15    lsm   created file
******************************************************************************/
#ifndef PARA_PADDS_ALGORITHM_H
#define PARA_PADDS_ALGORITHM_H

#include "MyHeaderInc.h"

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;

/******************************************************************************
class ParaPADDS
******************************************************************************/
class ParaPADDS : public AlgorithmABC
{
   public:
      ParaPADDS(ModelABC * pModel);
      ~ParaPADDS(void){ DBG_PRINT("ParaPADDS::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void){ return;}
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      int UpdateArchive(double * pX, int nX, double * pF, int nF);
      void DestroyArchive(ArchiveStruct * pArch);
      void F(ArchiveStruct * pA);
      void Calc_Z(ArchiveStruct * archive);
      void SortArchive(ArchiveStruct ** pArch, int size, int whichObj);
      void SortPoints(double ** X, int size, int which);
      int dominion_status(ArchiveStruct * x1, ArchiveStruct * x2);
      ArchiveStruct * SelectFrom(ArchiveStruct * pArchive);
      double neigh_value_continuous(double s, double s_min, double s_max, double r);
      double HV(int data_n, int dim_n, double * ref, double ** points);
      bool covers(double* cub, double * regLow);
      bool partCovers(double* cub, double * regUp);
      int containsBoundary(double* cub, double * regLow, int split);
      double getMeasure(double * regLow, double * regUp);
      int isPile(double* cub, double * regLow, double * regUp);
      double getMedian(double * bounds, int size);
      double computeTrellis(double * regLow, double * regUp, double * trellis);
      void stream(double * regionLow, double * regionUp, double ** points, int npoints, int split, double cover);
      int bool_vec_to_ulong(bool * pB, int size);
      void ulong_to_bool_vec(int val, bool * pB, int size);

      ModelABC * m_pModel;
      ArchiveStruct * m_pNonDom; //non-dominated solutions
      ArchiveStruct * m_pDom; //dominated solutions
      int m_NumNonDom; //number of non-dominated solutions
      int m_NumDom; //number of dominated solutions
      int m_CurIter;
      int m_nprocessors;
      int m_rank;
      int m_NumInit; //number of user-specified initial guesses
      double ** m_pInit; //list of user-specified initial guesses
      double * m_stest_flat; //params and objectives

      /* ==================================================================
      Store global variables from C++ version as private data members.
      ================================================================== */
      int m_maxiter; //number of total objective function evaluation
      int m_num_dec; //number of decesion variables (depends on the objective functions and the problem we are solving)
      int m_num_objs; //number of objective functions
      //selection metric
      //0: Random
      //1: Crowding distance
      //2: Hypervolume Contribution (ESTIMATE)
      //3: Hypervolume Contribution (EXACT)
      int m_Select_metric; 
      double m_fraction1;
      int m_dominance_flag;
      unsigned int m_seed;

      //member funcs for interfacing with hypervolume code
      int m_dim;
      int m_dimension;
      double m_dSqrtDataNumber;
      double m_volume;    
}; /* end class ParaPADDS */

extern "C" {
void PARA_PADDS_Program(int argC, StringType argV[]);
}

#endif /* PARA_PADDS_ALGORITHM_H */
