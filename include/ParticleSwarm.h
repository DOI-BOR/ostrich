/**************************************************************************************************************************************************************

Particle Swarm Optimization (PSO) applies concepts of social behavior to solve optimization problems. The PSO algorithm starts with a 'swarm' of particles
(solutions) and "flies" this population through the design space in search of the optimal solution. At each iteration, a given particle uses it's own prior
best solution (cognitive behavior) along with the current best solution of all particles (social behavior) to decide where to go next.

**************************************************************************************************************************************************************/

#ifndef PARTICLE_SWARM_H
#define PARTICLE_SWARM_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Include C++ classes
#include <algorithm>
#include <vector>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "QuadTree.h"
#include "StatUtility.h"
#include "LatinHypercube.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class ParticleSwarm : public Algorithm {
   public:
        // Define constructors and destructors
        ParticleSwarm();
        ~ParticleSwarm(void){ DBG_PRINT("ParticleSwarm::DTOR"); Destroy(); }
        void Destroy(void);

        // Define sampling functions used by the PSO algorithm
        void WarmStart(void);                                                    // todo: implement - placeholder for warm starting the swarm from a previous run
        std::vector<std::vector<double>> CreateInitialSwarm(int swarmSize);      // Function to create an initial swarm at the start of the analysis
        void UpdateSwarm(int generation);                                        // Function to update particle velocities/positions for the next iteration

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using PSO
        void Calibrate(void);                                               // Function to solve least squares fitting using PSO

        // Output functions
        void WriteStartingMetrics(void);                                      // Function to write the starting PSO metrics
        void WriteEndingMetrics(void);                                        // Function to write the ending PSO metrics


   private:
        // Configuration variables
        int m_SwarmSize = 20;                                                // Number of particles within the swarm
        int m_NumGenerationsMaximum = 50;                                    // Number of solution generations
        double m_Constriction = 1.00;                                        // Constriction factor
        double m_Cognitive = 2.00;                                           // Cognitive weight (c1)
        double m_Social = 2.00;                                              // Social weight (c2)
        double m_Inertia = 1.20;                                             // Inertia weight
        double m_InertiaInitial = 1.20;                                      // Starting inertia weight, retained so it can be restored/reduced each generation
        double m_InertiaReductionRate = 0.10;                                // Rate at which the inertia weight is reduced each generation
        bool m_bLinearInertiaReduction = false;                              // true: linearly reduce inertia weight to zero over the course of the run
        int m_InitType = RANDOM_INIT;                                        // Initializaiton strategy used to create the initial swarm
        double m_StopVal = 0.001;                                            // Convergence tolerance used to stop the analysis

        // Solution tracking variables
        int m_Generation = 0;                                                // Counter for the current generation number
        double m_CurStop = 1000.0;                                           // Current convergence value to be compared against m_StopVal

        double m_BestObjectiveIteration = INFINITY;                          // Best objective within the iteration
        int m_BestObjectiveIndexIteration = -1;                              // Best objective index within the iteration

        int m_NumUprViols = 0;                                               // Counter for upper bound violations (used in reporting)
        int m_NumLwrViols = 0;                                               // Counter for lower bound violations (used in reporting)

        // Quad-tree initialization support
        QuadTree* m_pTrees = NULL;
        int m_TreeSize = 0;

        // Swarm state - one entry per particle
        std::vector<std::vector<double>> m_Positions;                        // Current position of each particle
        std::vector<std::vector<double>> m_Velocities;                       // Current velocity of each particle
        std::vector<std::vector<double>> m_PersonalBestPositions;            // Personal (local) best position of each particle
        std::vector<double> m_PersonalBestObjectives;                        // Personal (local) best objective of each particle

        // Solution functions
        double CalcSwarmMedian(std::vector<double> objectives);              // Function to calculate the median objective of the swarm
        void ConstrainToBounds(int particleIndex);                           // Function to constrain a particle to stay within parameter bounds, preserving direction of movement

}; 

extern "C" {
void PSO_Program(int argC, StringType argV[]);
}

#endif 
