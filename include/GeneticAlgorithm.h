/******************************************************************************
File     : GeneticAlgorithm.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

A Genetic Algorithm applies concepts (namely survival of the fittest and 
natural selection) from evolutionary theory to optimization problems. The 
Genetic Algorithm starts with a population of coded solutions (ChromosomePool) 
and evolves this population using the processes of Selection, Crossover and 
Mutation such that each successive genration of solutions is an improvement 
(on average) over previous generations.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
******************************************************************************/

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class StatsClass;
class ModelABC;
class ChromosomePool;

/******************************************************************************
class GeneticAlgorithm

******************************************************************************/
class GeneticAlgorithm : public AlgorithmABC
{
   public:
      GeneticAlgorithm(ModelABC * pModel);
      ~GeneticAlgorithm(void){ DBG_PRINT("GeneticAlgorithm::DTOR"); Destroy(); }
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurGen; }

   private:
      ModelABC * m_pModel;
      ChromosomePool * m_pPopulation;
      StatsClass * m_pStats;
      double m_StopVal; //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)
      int m_MaxGens;
      int m_CurGen;
}; /* end class GeneticAlgorithm */

extern "C" {
void GA_Program(int argC, StringType argV[]);
}

#endif /* GENETIC_ALGORITHM_H */


