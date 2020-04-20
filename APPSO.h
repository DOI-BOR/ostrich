/******************************************************************************
File     : APPSO.h
Author   : L. Shawn Matott
Copyright: 2014, L. Shawn Matott

Asynchronous Parallel Particle Swarm Optimization (APPSO).

A parallel version of PSO based on the asynchronous master-slave approach
of the PDDS algorithm.

Version History
12-30-14    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef APPSO_H
#define APPSO_H

#include "MyHeaderInc.h"

//paraent class
#include "AlgorithmABC.h"

//forward declarations
class ModelABC;
class StatsClass;

/******************************************************************************
class APPSO

******************************************************************************/
class APPSO : public AlgorithmABC
{
   public:
      APPSO(ModelABC * pModel);
      ~APPSO(void){ DBG_PRINT("APPSO::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurGen; }

   private:
      void EvaluateSwarm(int id, int nprocs, int gen);
      double CalcPSOMedian(void);
      void MakeParameterCorrections(double * x, double * xb, int n, double a);

      ModelABC * m_pModel;
      ParticleStruct * m_pSwarm;
      StatsClass * m_pStats;
      int m_SwarmSize;
      int m_MaxGens;
      int m_BestIdx;
      double m_Best;
      double m_Constrict; //constriction factor
      double m_c1; //cognitive weight
      double m_c2; //social weight
      double m_Inertia; //Inertia weight
      double m_RedRate;   
      bool m_LinRedFlag; //true: linearly reduce interia to zero
      int m_CurGen;
      int m_id;
      int * m_Assignments; //keeps track of which particle each processor is currently working on    

      //buffer for initial parameter values
      int       m_NumInit;
      double ** m_pInit;

      //metrics
      int m_NumUprViols;
      int m_NumLwrViols;
      double * m_Fmedian;
}; /* end class ParticleSwarm */

extern "C" {
void APPSO_Program(int argC, StringType argV[]);
}

#endif /* APPSO_H */


