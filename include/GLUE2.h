/**************************************************************************************************************************************************************

Generalized Likelihood Uncertainty Estimation (GLUE), as described in Beven and Binley (1992). GLUE repeatedly draws batches of uniformly random parameter
sets and retains the best NumBehavioral samples. Samples with an objective below the behavioral threshold are considered behavioral. The search ends when
the desired number of behavioral samples has been found, or when the sample budget (MaxSamples) has been consumed. Each batch of samples is solved using
the MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

#ifndef GLUE_H
#define GLUE_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

// Include C++ classes
#include <vector>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "ParameterGroup.h"
#include "ReadUtility.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class GLUE : public Algorithm {
    public:
        // Define constructors and destructors
        GLUE();
        ~GLUE(void) { DBG_PRINT("GLUE::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to approximate uncertainty using GLUE

        // Output functions
        void WriteStartingMetrics(void);                                    // Function to write the starting GLUE metrics
        void WriteEndingMetrics(void);                                      // Function to write the ending GLUE metrics

    private:
        // A retained sample
        struct Sample {
            double objective = HUGE_VAL;                                    // Objective function value
            std::vector<double> parameters;                                 // Parameter set (transformed). Empty until a sample is retained.
        };

        // Configuration variables
        int m_SamplesPerIter = 10;                                          // Number of samples generated and solved per iteration
        int m_NumDesired = 10;                                              // Desired number of behavioral samples
        long long m_MaxSamples = 100;                                       // Sample budget
        double m_Threshold = 1000.00;                                       // Behavioral threshold (behavioral samples are strictly below)

        // Solution tracking variables
        std::vector<Sample> m_Behavioral;                                   // The best m_NumDesired samples, sorted from best to worst
        int m_NumFound = 0;                                                 // Number of retained samples that are behavioral
        int m_MaxIters = 0;                                                 // Maximum number of iterations
        int m_CurIter = 0;                                                  // Number of iterations performed

        // Solution functions
        std::vector<double> GenerateSample(void);                           // Function to generate a uniformly random parameter set
        void UpdateBehavioral(std::vector<std::vector<double>>& samples, std::vector<double>& objectives); // Function to revise the retained samples
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group
};

extern "C" {
void GLUE_Program(int argC, StringType argV[]);
}

#endif /* GLUE_H */
