/**************************************************************************************************************************************************************

DDS for Approximation of Uncertainty (DDS-AU), as described in Tolson and Shoemaker (2008), Water Resources Research. DDS-AU performs a number of independent,
randomly initialized, short-budget DDS searches. A single parameter set is selected from each search (either the best solution, or a randomly selected
behavioral solution), and the selected parameter sets with an objective at or below the behavioral threshold approximate the parameter uncertainty.

The pre-refactor implementation ran each search with a separate DDS (or parallel DDS) algorithm object. The Algorithm base class configures and terminates the
MPI workers once per analysis, so this version performs the DDS searches within a single algorithm. Each search follows the DDS algorithm exactly. When
ParallelSearches is enabled, the searches are advanced together: one candidate from every active search is solved concurrently at each step. Otherwise, the
searches are performed one after another, and only the random initialization samples of each search are solved concurrently.

**************************************************************************************************************************************************************/

#ifndef DDSAU_ALGORITHM_H
#define DDSAU_ALGORITHM_H

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


class DDSAU : public Algorithm {
    public:
        // Define constructors and destructors
        DDSAU();
        ~DDSAU(void) { DBG_PRINT("DDSAU::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to approximate uncertainty using DDS-AU

        // Output functions
        void WriteStartingMetrics(void);                                    // Function to write the starting DDS-AU metrics
        void WriteEndingMetrics(void);                                      // Function to write the ending DDS-AU metrics

    private:
        // Information tracked for each independent DDS search
        struct Search {
            int budget = 0;                                                 // Number of objective function evaluations allocated to the search
            int initFunctEvals = 0;                                         // Number of random initialization samples
            bool bRevised = false;                                          // true: results were read from a previous analysis (ReviseAU)
            double Fbest = INFINITY;                                        // Objective of the current best solution
            std::vector<double> bestParams;                                 // Current best solution (transformed)
            std::vector<double> objectives;                                 // Objective of every evaluation, in order
            std::vector<std::vector<double>> parameters;                    // Parameter set of every evaluation (transformed), in order
            std::vector<int> recordEvals;                                   // Evaluation number of each improvement to the best solution
            int selected = -1;                                              // Index of the evaluation selected to represent the search
        };

        // Configuration variables
        double m_r_val = 0.2;                                               // Perturbation value (0 < r < 1), as a fraction of the parameter range
        int m_nsols = 25;                                                   // Number of independent searches
        int m_MinIter = 30;                                                 // Minimum number of evaluations per search
        int m_MaxIter = 70;                                                 // Maximum number of evaluations per search
        bool m_bParallel = false;                                           // true: advance the searches together
        double m_fmax = 1000.00;                                            // Behavioral threshold
        bool m_bRandomize = false;                                          // true: randomly select a behavioral solution from each search
        bool m_bReviseAU = false;                                           // true: re-use the results of previously completed searches

        // Solution tracking variables
        std::vector<Search> m_Searches;                                     // Independent DDS searches
        int m_nbhvr = 0;                                                    // Number of searches whose selected solution is behavioral
        int m_CurIter = 0;                                                  // Number of objective function evaluations performed

        // Solution functions
        void RunSearches(std::vector<int> searchIndices);                   // Function to perform a set of DDS searches together
        void SelectSolution(Search& search);                                // Function to select the solution that represents a search
        bool ReadSearch(int searchIndex);                                   // Function to read the results of a previous search (ReviseAU)
        std::vector<double> Evaluate(std::vector<std::vector<double>> samples); // Function to solve a set of parameter sets and track the best solution
        std::vector<double> CorrectParameters(std::vector<double> parameters); // Function to apply any parameter corrections
        double PerturbParam(double x_best, ParameterABC* pParam);           // Function to generate a neighboring value for a single parameter
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group

        // Output functions
        void WriteSearch(int searchIndex);                                  // Function to write the output files of a search
        void WriteSummary(void);                                            // Function to write the OstOutputDDSAU.txt summary file
        void WriteSolutionLine(FILE* pFile, int run, double objective, std::vector<double> parameters); // Function to write an evaluation
};

extern "C" {
void DDSAU_Program(int argC, StringType argV[]);
}

#endif /* DDSAU_ALGORITHM_H */
