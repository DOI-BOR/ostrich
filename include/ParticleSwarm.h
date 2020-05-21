/******************************************************************************
File     : ParticleSwarm.h
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Particle Swarm Optimization (PSO) applies concepts of social behavior to
solve optimization problems. The PSO algorithm starts with a 'swarm' of 
particles (solutions) and "flies" this population through the design space in
search of the optimal solution. At each iteration, a given particle uses it's
own prior best solution (cognitive behavior) along with the current best 
solution of all particles (social behavior) to decide where to go next.

Version History
02-25-04    lsm   added copyright information and initial comments.
03-24-04    lsm   added hybrid of PSO and levenberg-marquardt
07-08-04    lsm   added parallel support
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
11-18-04    lsm   Added convergence criteria, based on median fitness of swarm.
12-02-04    lsm   Added support for the seeding of particle swarm
04-14-05    lsm   Added support for linearly reducing the inertia weight to zero.
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
******************************************************************************/

#ifndef PARTICLE_SWARM_H
#define PARTICLE_SWARM_H

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;
class StatsClass;
class QuadTree;

/******************************************************************************
class ParticleSwarm

******************************************************************************/
class ParticleSwarm : public AlgorithmABC
{
   public:
      ParticleSwarm(ModelABC * pModel);
      ~ParticleSwarm(void){ DBG_PRINT("ParticleSwarm::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurGen; }

   private:
      void EvaluateSwarm(void);
      void BcastSwarm(void);
      void EvalSwarmParallel(void);
      void EvalSwarmSuperMUSE(void);
      double CalcPSOMedian(void);
     void MakeParameterCorrections(double * x, double * xb, int n, double a);

      ModelABC * m_pModel;
      ParticleStruct * m_pSwarm;
      StatsClass * m_pStats;
      QuadTree * m_pTrees;
      int m_TreeSize;
      int m_SwarmSize;
      int m_MaxGens;
      int m_BestIdx;
      double m_Best;
      double m_Constrict; //constriction factor
      double m_c1; //cognitive weight
      double m_c2; //social weight
      double m_Inertia; //Inertia weight
      double m_RedRate;   
      PopInitType m_InitType;      
      bool m_LinRedFlag; //true: linearly reduce interia to zero
      int m_CurGen;
      double m_StopVal;  //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)

      //buffers used in MPI-parallel communication
      double * m_pBuf;
      double * m_pMyBuf;
      double * m_pTmpBuf;
      double * m_pBigBuf;

      //buffer for initial parameter values
      int       m_NumInit;
      double ** m_pInit;

      //metrics
      int m_NumUprViols;
      int m_NumLwrViols;
      double * m_Fmedian;
}; /* end class ParticleSwarm */

extern "C" {
void PSO_Program(int argC, StringType argV[]);
void PSO_LEVMAR_Program(int argC, StringType argV[]);
}

#endif /* PARTICLE_SWARM_H */


