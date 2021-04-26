/******************************************************************************

A Genetic Algorithm applies concepts (namely survival of the fittest and 
natural selection) from evolutionary theory to optimization problems. The 
Genetic Algorithm starts with a population of coded solutions (ChromosomePool) 
and evolves this population using the processes of Selection, Crossover and 
Mutation such that each successive genration of solutions is an improvement 
(on average) over previous generations.

******************************************************************************/

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include "MyHeaderInc.h"

//parent class
#include "Algorithm.h"
#include "Gene.h"
#include "QuadTree.h"
#include "StatUtility.h"
#include "ModelWorker.h"
#include "WriteUtility2.h"


//forward declarations
//class StatsClass;

/******************************************************************************
class GeneticAlgorithm

******************************************************************************/
class GeneticAlgorithm : Algorithm {
   public:
      GeneticAlgorithm();      
      ~GeneticAlgorithm(void){ DBG_PRINT("GeneticAlgorithm::DTOR"); Destroy(); }
      double** CreateInitialSample(int sampleSize);
      double** CreateSample(double* objectives, int numberOfObjectives, double** samples);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteStartingMetrics(void);
      void WriteEndingMetrics(void);
      void WarmStart(void);

   private:
       // Configuration variables
       int m_NumGenerationsMaximum = 10;            // Number of solution generations 
       int m_NumSurvivors = 1;                      // Number of best unmodified survivors per generation
       int m_NumPopulation = 50;                    // Number of members within each generation
       double m_MutationRate = 0.05;                // Mutation rate within the population       
       double m_CrossoverRate = 0.90;               // Crossover rate within the population
       int m_InitType = RANDOM_INIT;                // Initializaiton strategy used to create the initial generation 
       double m_StopVal = 0.0001;                   // Convergence tolerance used to stop the analysis
       double m_CurStop = 1000.0;                   //current convergence val (compared against m_StopVal)
                                                    // Solution variables
       int m_Generation = 0;                        // Counter for the current generation number
       RealEncodedGene *m_pGenes;                   // Genes that represent the parameters

       double m_BestObjectiveIteration = INFINITY;       // Best objective within the iteration
       int m_BestObjectiveIndexIteration = -1;      // Best objective index within the iteration
       
       void GetBestObjective(double* objectives, int numberOfObjectives);
       double CalcMeanFitness(double* objectives, int numberOfObjectives);
       double CalcMedianFitness(double* objectives, int numberOfObjectives);

       void TourneySelection(int nCombatants, double* objectives, int numberOfObjectives, double** samples, double** scratch);
       void Crossover(double* objectives, int numberOfObjectives, double** samples, double** scratch);
       void Mutate(double** scratch, int numberOfSamples);

       //StatsClass m_pStats;

      
}; 

extern "C" {
void GA_Program(int argC, StringType argV[]);
}

#endif /* GENETIC_ALGORITHM_H */


