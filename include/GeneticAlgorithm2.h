/**************************************************************************************************************************************************************

A Genetic Algorithm applies concepts (namely survival of the fittest and natural selection) from evolutionary theory to optimization problems. The
Genetic Algorithm starts with a population of coded solutions and evolves this population using the processes of Crossover and  Mutation such that each
successive genration of solutions is an improvement  (on average) over previous generations.

**************************************************************************************************************************************************************/

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Include C++ classes
#include <algorithm>
#include <numeric>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "Gene.h"
#include "QuadTree.h"
#include "StatUtility.h"
#include "LatinHypercube.h"
#include "ParameterGroup.h"
#include <ParameterABC.h>
//#include "StatsClass.h"
#include "Exception.h"
#include "Utility.h"



class GeneticAlgorithm : public Algorithm {
   public:
        // Define constructors and destructors
        GeneticAlgorithm();      
        ~GeneticAlgorithm(void){ DBG_PRINT("GeneticAlgorithm::DTOR"); Destroy(); }
        void Destroy(void);

        // Define sampling functions used by the GA algorithm
        void WarmStart(void);                                                   // Function to start from a previous solution
        std::vector<std::vector<double>> CreateInitialSample(int sampleSize);   // Function to create an initial sample at the start of the analysis
        void CreateSample(std::vector<double>& objectives, std::vector<double>& objectivesScratch, std::vector<std::vector<double>>& samples, 
                        std::vector<std::vector<double>>& samplesScratch);      // Function to create ongoing samples for each iteration

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using GA
        void Calibrate(void);                                               // Function to solve least squares fitting using GA
      
        // Output functions
        void WriteStartingMetrics(void);                                      // Function to write the starting GA metrics
        void WriteEndingMetrics(void);                                        // Function to write the ending GA metrics
      

   private:
        // Configuration variables
        int m_NumGenerationsMaximum = 10;                                    // Number of solution generations 
        int m_NumSurvivors = 1;                                              // Number of best unmodified survivors per generation
        int m_NumPopulation = 50;                                            // Number of members within each generation
        double m_MutationRate = 0.05;                                        // Mutation rate within the population       
        double m_CrossoverRate = 0.90;                                       // Crossover rate within the population
        int m_InitType = RANDOM_INIT;                                        // Initializaiton strategy used to create the initial generation 
        double m_StopVal = 0.0001;                                           // Convergence tolerance used to stop the analysis
        RealEncodedGene* m_pGenes;                                           // Genes that represent the parameters

        // Solution tracking variables
        int m_Generation = 0;                                                // Counter for the current generation number
        double m_CurStop = 1000.0;                                           // Current convergence value to be compared against m_StopVal
        
        double m_BestObjectiveIteration = INFINITY;                          // Best objective within the iteration
        int m_BestObjectiveIndexIteration = -1;                              // Best objective index within the iteration
       
        // Solution functions
        void GetBestObjective(std::vector<double> objectives);              // Function to get the best objective from the population
        double CalcMeanFitness(std::vector<double> objectives);             // Function to calculate the mean fitness of the population
        double CalcMedianFitness(std::vector<double> objectives);           // Function to calcualte the median fitness of the population

        void Crossover(std::vector<double>& objectives, std::vector<double>& objectivesScratch, std::vector<std::vector<double>>& samples,
                      std::vector<std::vector<double>>& samplesScratch);    // Function to crossover the population to generate a new population
        void Mutate(std::vector<std::vector<double>>& scratch);             // Function to mutate the population to generate a new population
      
}; 

extern "C" {
void GA_Program(int argC, StringType argV[]);
}

#endif 


