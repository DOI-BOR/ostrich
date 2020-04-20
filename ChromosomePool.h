/******************************************************************************
File     : ChromosomePool.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The ChromosomePool is a container for a set of Chromosomes (coded design 
variables). The Genetic Algorithm uses ChromosomePools to store the population 
of solutions for a given generation along with the mating pool from which the 
next generation of solutions will be produced.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes
                  Added metrics collection and reporting
                  Added intermediate output (when generation is being evaluated)
11-18-04    lsm   Added convergence criteria, based on median fitness of population.
03-21-05    lsm   Added support for seeding initial population
01-01-07    lsm   Algorithms now use an abstract model base class (ModelABC).
07-18-07    lsm   Added SuperMUSE support
01-13-15    lsm   Added support for asynchrounous parallel
******************************************************************************/
#ifndef CHROMOSOME_POOL_H
#define CHROMOSOME_POOL_H

#include "MyHeaderInc.h"

//forward declarations
class ChromosomeCommunicator;
class Chromosome;
class ModelABC;
class QuadTree;
class ParameterGroup;

#define PARALLEL_TYPE_SYNCH  (0)
#define PARALLEL_TYPE_ASYNCH (1)

/******************************************************************************
class ChromosomePool

 class models the collection of chromosome.
******************************************************************************/
class ChromosomePool
{
   public:
      ChromosomePool(void);
      ~ChromosomePool(void){ DBG_PRINT("ChromosomePool::DTOR"); Destroy();}
      void Destroy(void);

      int GetPoolSize(void){ return m_PoolSize;}
      void EvalFitness(void);
      void CreateNxtGen(void);
      void CreateNxtGen(double pct);
      double CalcAvgFitness(void);
      double CalcMedianFitness(void);
      Chromosome * GetBestFit(void);
      void Initialize(void);
      void Initialize(int * budget);
      void CreateComm(ModelABC * pModel);      
      ParameterGroup * ConvertChromosome(Chromosome * pChromo);
      int GetNumGens(void){return m_NumGenerations;}
      double GetStopVal(void){return m_StopVal;}
      void WriteMetrics(FILE * pFile);
      void SetChromosome(int i, double * vals);

   private:
      void TourneySelection(int nCombatants);
      void Crossover(void);
      void Mutate(void);
      void CreatePrototype(void);
      void EvalFitParallel(void);
      void EvalFitSuperMUSE(void);
      void BcastPopulation(void);
      void FreezeGenes(int numFreeze);        
      void EvalFitnessAsynch(int rank, int nprocs);

      QuadTree * m_pTrees;
      int m_TreeSize;
      /*--------------------------------------
      Two pools, one for current generation 
      and a scratch pool for creating the next 
      generation
      ---------------------------------------*/
      Chromosome **  m_pPool;
      Chromosome **  m_pScratch;
      int m_PoolSize;
      int m_NumInit;
      double ** m_pInit;

      ChromosomeCommunicator * m_pComm;
      Chromosome * m_Proto; 
      int m_Generation;
      int m_NumSurvivors;
      int m_NumGenerations;
      PopInitType m_InitType;

      //type of parallization
      // 0 = synchronous
      // 1 = asynchronous
      int m_ParallelType;

      //read from file by the pool, but passed up to the GeneticAlg parent class
      double m_StopVal;

      //buffers used in MPI-parallel communication
      double * m_pBuf;      
      double * m_pMyBuf;
      double * m_pTmpBuf;
      double * m_pBigBuf;

      //list of chromosomes currently assigned to each slave
      int * m_pAssignments;

      //metrics
      int * m_pMutCount;
      double * m_Fmedian;
}; /* end class ChromosomePool */

#endif /* CHROMOSOME_POOL_H */



