/**************************************************************************************************************************************************************

Balanced Exploration-Exploitation Random Search (BEERS) applies a balanced approach to randomly search a parameter space. The search initially favors
exploration of the search space and will try to maximize the distance between points that are evaluated. As the search progresses the algorithm favors
exploitation and will favor points that are close to the current optimal. To facilitate exploration vs. exploitation the algorithm maintains an archive of
every point that is ever evaluated during the search.

Each candidate depends on the objective of every previously evaluated candidate, so the samples are evaluated one at a time using the shared MPI
work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

#ifndef BEERS_H
#define BEERS_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

// Include C++ classes
#include <vector>
#include <string>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "ParameterGroup.h"
#include "ReadUtility.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class BEERS : public Algorithm {
   public:
        // Define constructors and destructors
        BEERS();
        ~BEERS(void){ DBG_PRINT("BEERS::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to perform the BEERS search
        void Calibrate(void);                                               // Function to perform the BEERS search in LS fitting mode

        // Output functions
        void WriteStartingMetrics(void);                                     // Function to write the starting BEERS metrics
        void WriteEndingMetrics(void);                                       // Function to write the ending BEERS metrics

   private:
        // Archive member
        struct ArchiveMember {
            std::vector<double> X;                                          // Parameters (transformed)
            double F = INFINITY;                                            // Objective
            double Z = 0;                                                   // Rank: number of archive members that are better or equal
            double P = 0;                                                   // Model probability
        };

        // Configuration variables
        int m_NumSamples = 25;                                              // Number of samples to evaluate (computational budget)

        // Search space information
        std::vector<double> m_Min;                                          // Lower bound (transformed) for each parameter
        std::vector<double> m_Range;                                        // Range (transformed) for each parameter
        double m_MaxDist = 1.00;                                            // Maximum normalized distance in the search space
        double m_MinProbAccept = 0.00;                                      // Minimum acceptance probability, increased to guarantee acceptance

        // Solution tracking variables
        std::vector<ArchiveMember> m_Archive;                               // Every evaluated sample, in evaluation order
        int m_BestIndex = 0;                                                // Index of the best archive member
        int m_CurSample = 0;                                                // Number of samples evaluated

        // Solution functions
        double Evaluate(std::vector<double> parameters);                    // Evaluates a single sample
        std::vector<double> ApplyParameterGroup(std::vector<double> parameters); // Sets values into the parameter group and reads them back (integer snapping)
        void SetParameterGroup(std::vector<double> parameters);             // Sets values into the parameter group
        void AssignModelProbs(double Fbest);                                // Computes the likelihood of each archive member
        void CalcProbabilities(const std::vector<double>& candidate, double& Pexploit, double& Pexplore);
        void AdjustRanks(int currentIndex);                                 // Adjusts the archive ranks for the newest member
        void WriteArchive(void);                                            // Writes the archive to OstArchive.txt
};

extern "C" {
void BEERS_Program(int argC, StringType argV[]);
}

#endif /* BEERS_H */
