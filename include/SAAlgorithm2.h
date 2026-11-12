/**************************************************************************************************************************************************************

Simulated Annealing (SA) for continuous parameters. The design space is first "melted" with a random walk to estimate the initial temperature. The system is
then cooled over a series of outer iterations. Within each outer iteration, a number of transitions (moves to neighboring solutions) are attempted and are
accepted if they improve the objective or, otherwise, with a probability that decreases as the temperature drops (Metropolis criterion). Transitions are
generated in batches that are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

#ifndef SA_ALGORITHM_H
#define SA_ALGORITHM_H

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


class SAAlgorithm : public Algorithm {
   public:
        // Define constructors and destructors
        SAAlgorithm();
        ~SAAlgorithm(void){ DBG_PRINT("SAAlgorithm::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using SA

        // Output functions
        void WriteStartingMetrics(void);                                     // Function to write the starting SA metrics
        void WriteEndingMetrics(void);                                       // Function to write the ending SA metrics

   private:
        // Configuration variables
        int m_MaxOuter = 20;                                                // Maximum number of outer iterations (temperature reductions)
        int m_MaxInner = 10;                                                // Number of transitions attempted per outer iteration
        int m_NumMelts = 100;                                               // Number of evaluations used to determine the initial temperature
        double m_TempFactor = 0.90;                                         // Temperature reduction factor
        double m_FinalTemp = 10.00;                                         // Final temperature (user-specified or computed)
        TempMethodType m_FinalTempMethod = TMETHD_NORM;                     // How the final temperature is determined
        TransMethodType m_TransitionMethod = TRANS_GAUSS;                   // How transitions are generated
        double m_StopVal = 0.001;                                           // Convergence tolerance used to stop the analysis

        // Solution tracking variables
        int m_NumOuter = 0;                                                 // Counter for the current outer iteration number
        double m_InitTemp = 10.00;                                          // Initial temperature, computed by the melting operation
        double m_CurTemp = 10.00;                                           // Current temperature
        double m_dEavg = 0.00;                                              // Average energy change, estimated from the melting operation
        double m_CurStop = 1.00;                                            // Current convergence value to be compared against m_StopVal
        std::vector<double> m_Current;                                      // Current (most recently accepted) parameter set
        int m_BatchSize = 1;                                                // Number of transitions evaluated concurrently

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
        std::vector<double> UniformTransition(void);                        // Function to generate a uniform transition from the current point
        std::vector<double> GaussTransition(double currentObjective);       // Function to generate a Gaussian transition from the current point
        double AcceptTransition(std::vector<double> proposal, double proposalObjective, double currentObjective); // Metropolis acceptance test
        void GenerateRandomMove(std::vector<double>& parameters);           // Function to randomly displace a parameter set (melting)
        std::vector<double> CorrectParameters(std::vector<double> parameters); // Function to apply any parameter corrections
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group

};

extern "C" {
void SA_Program(int argC, StringType argV[]);
}

#endif
