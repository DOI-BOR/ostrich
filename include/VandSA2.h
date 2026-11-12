/**************************************************************************************************************************************************************

Simulated Annealing (SA) based on: Vanderbilt and Louie. 1984. "A Monte Carlo Simulated Annealing Approach to Optimization over Continuous Variables".
Journal of Computational Physics. vol. 56, pg. 259-271. The design space is first "melted" with a random walk to estimate the initial temperature. The system
is then cooled over a series of outer iterations. Within each outer iteration, a number of transitions are attempted and accepted according to the Metropolis
criterion. With the Vanderbilt-Louie transition method, the step distribution is adapted after each outer iteration using the covariance of the accepted
parameter sets. Transitions are generated in batches that are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm
base class.

**************************************************************************************************************************************************************/

#ifndef VAND_SA_H
#define VAND_SA_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

// Include C++ classes
#include <algorithm>
#include <vector>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "StatUtility.h"
#include "ParameterGroup.h"
#include "ReadUtility.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class VandSA : public Algorithm {
   public:
        // Define constructors and destructors
        VandSA();
        ~VandSA(void){ DBG_PRINT("VandSA::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using Vanderbilt-Louie SA

        // Output functions
        void WriteStartingMetrics(void);                                     // Function to write the starting SA metrics
        void WriteEndingMetrics(void);                                       // Function to write the ending SA metrics

   private:
        // Configuration variables
        int m_MaxOuter = 20;                                                // Maximum number of outer iterations (temperature reductions)
        int m_MaxInner = 10;                                                // Number of transitions attempted per outer iteration
        int m_NumMelts = 100;                                               // Number of evaluations used to determine the initial temperature
        double m_TempFactor = 0.90;                                         // Temperature reduction factor
        double m_FinalTemp = 10.00;                                         // User-specified final temperature
        bool m_bUserFinalTemp = false;                                      // true: the user has specified the desired final temperature
        TransMethodType m_TransitionMethod = TRANS_GAUSS;                   // How transitions are generated
        double m_StopVal = 0.001;                                           // Convergence tolerance used to stop the analysis

        // Solution tracking variables
        int m_NumOuter = 0;                                                 // Counter for the current outer iteration number
        double m_InitTemp = 10.00;                                          // Initial temperature, computed by the melting operation
        double m_CurTemp = 10.00;                                           // Current temperature
        double m_CurStop = 1.00;                                            // Current convergence value to be compared against m_StopVal
        std::vector<double> m_Current;                                      // Current (most recently accepted) parameter set
        int m_BatchSize = 1;                                                // Number of transitions evaluated concurrently

        // Matrices needed for the implementation of Vanderbilt and Louie (1984)
        std::vector<std::vector<double>> m_Q;                               // Step distribution matrix
        std::vector<std::vector<double>> m_QT;                              // Q transpose
        std::vector<std::vector<double>> m_Cov;                             // Covariance matrix ('s')
        std::vector<std::vector<double>> m_Shape;                           // Random walk shape matrix ('S')

        // Metrics
        int m_MeltCount = 0;                                                // Number of melting evaluations
        int m_TransCount = 0;                                               // Number of transition evaluations
        int m_EquilCount = 0;                                               // Number of equilibrations
        int m_NumAborts = 0;                                                // Number of rejected transitions
        int m_NumUprViols = 0;                                              // Counter for upper bound violations
        int m_NumLwrViols = 0;                                              // Counter for lower bound violations
        int m_NumUphill = 0;                                                // Number of accepted uphill transitions
        int m_NumDownhill = 0;                                              // Number of accepted downhill transitions
        double m_CurProb = -1.00;                                           // Average acceptance probability of the latest equilibration
        double m_InitProb = -1.00;                                          // Average acceptance probability of the first equilibration
        double m_TotProb = 0.00;                                            // Accumulated acceptance probability within an equilibration
        int m_NumProbTests = 0;                                             // Number of acceptance probability tests within an equilibration

        // Solution functions
        double Initialize(void);                                            // Function to evaluate the starting point and nearby perturbations
        double Melt(double bestObjective);                                  // Function to melt the design space and set the temperature schedule
        double Equilibrate(double bestObjective);                           // Function to perform the transitions at the current temperature
        std::vector<double> VanderbiltTransition(void);                     // Function to generate a Vanderbilt-Louie transition from the current point
        std::vector<double> GaussTransition(double currentObjective);       // Function to generate a Gaussian transition from the current point
        double AcceptTransition(std::vector<double> proposal, double proposalObjective, double currentObjective); // Metropolis acceptance test
        void UpdateStepMatrix(std::vector<std::vector<double>> transitions); // Function to update the step distribution matrix (Q)
        void GenerateRandomMove(std::vector<double>& parameters);           // Function to randomly displace a parameter set (melting)
        std::vector<double> CorrectParameters(std::vector<double> parameters); // Function to apply any parameter corrections
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group

};

extern "C" {
void VSA_Program(int argC, StringType argV[]);
}

#endif
