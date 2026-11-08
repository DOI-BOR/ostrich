/**************************************************************************************************************************************************************

The PADDS (Pareto Archived Dynamically Dimensioned Search) algorithm is a multi-objective version of the DDS algorithm (Asadzadeh and Tolson, 2009, 2013).
It was ported to OSTRICH by L. Shawn Matott from a C++ implementation of Mohammadamin Jahanpour (mjahanpo@uwaterloo.ca).

PADDS maintains an archive of non-dominated and dominated solutions. The initial solutions (user-specified and random) are solved together using the MPI
work-farming functionality provided by the Algorithm base class. Each subsequent candidate is generated from a member of the non-dominated archive and
depends on the result of the previous candidate, so the main search solves one candidate at a time. The objectives are the WSSE of each observation group.

**************************************************************************************************************************************************************/

#ifndef PADDS_ALGORITHM_H
#define PADDS_ALGORITHM_H

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


class PADDS : public Algorithm {
    public:
        // Define constructors and destructors
        PADDS();
        ~PADDS(void) { DBG_PRINT("PADDS::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to search for the Pareto front using PADDS

        // Output functions
        void WriteStartingMetrics(void);                                    // Function to write the starting PADDS metrics
        void WriteEndingMetrics(void);                                      // Function to write the ending PADDS metrics

    protected:
        // Constructor for derived algorithms (e.g. parallel PADDS) that read their own configuration
        PADDS(bool bReadConfiguration);

        // An archived solution
        struct Solution {
            std::vector<double> X;                                          // Parameter set (transformed)
            std::vector<double> F;                                          // Objectives
            double Z = -999.999;                                            // Selection metric
        };

        // Algorithm description
        std::string m_SetupName = "PADDS";                                  // Name written to the setup section of OstOutput0.txt
        std::string m_AlgorithmName = "PADDS - Pareto Archived Dynamically Dimensioned Search"; // Name written to the algorithm metrics
        int m_BatchSize = 0;                                                // Number of candidates solved together (0: not reported)

        // Configuration variables
        int m_maxiter = 50;                                                 // Number of objective function evaluations
        double m_fraction1 = 0.2;                                           // Perturbation value
        int m_Select_metric = 3;                                            // 0: Random, 1: Crowding distance, 2: Hypervolume contribution (estimate), 3: (exact)
        std::vector<std::vector<double>> m_Init;                            // User-specified initial parameter sets (transformed)

        // Solution tracking variables
        std::vector<Solution> m_NonDom;                                     // Non-dominated solutions. New solutions are inserted at the front.
        std::vector<Solution> m_Dom;                                        // Dominated solutions. New solutions are inserted at the front.
        int m_num_dec = 0;                                                  // Number of decision variables
        int m_num_objs = 0;                                                 // Number of objectives
        int m_dominance_flag = 0;                                           // -1: the last candidate was dominated by its parent
        int m_CurIter = 0;                                                  // Number of objective function evaluations performed
        int m_NumDuplicates = 0;                                            // Number of candidates discarded as duplicates
        unsigned int m_seed = 0;                                            // Random seed

        // Hypervolume variables
        int m_dimension = 0;
        double m_dSqrtDataNumber = 0.0;
        double m_volume = 0.0;

        // Bounds of the decision variables (transformed)
        std::vector<double> m_S_min;
        std::vector<double> m_S_max;

        // Configuration functions
        void ReadConfiguration(std::vector<std::string> blockNames, bool bReadInitParams); // Function to read the algorithm block and initial sets
        void BeginAnalysis(void);                                           // Function to validate the setup, configure the workers and write the setup

        // Search functions
        std::vector<double> RandomSolution(void);                           // Function to generate a uniformly random parameter set
        Solution SelectParent(void);                                        // Function to select the solution to perturb
        Solution PerturbSolution(Solution& sbest, double Pn);               // Function to generate a candidate from a selected solution
        bool ProcessCandidate(Solution& stest, Solution& sbest);            // Function to add a solved candidate to the archives

        // Solution functions
        std::vector<std::vector<double>> Evaluate(std::vector<std::vector<double>> samples); // Function to solve a set of parameter sets
        std::vector<double> ApplyParameterGroup(std::vector<double> parameters); // Function to set values into the parameter group and read them back
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group
        int UpdateArchive(std::vector<double> X, std::vector<double> F);    // Function to add a solution to the archives
        bool IsDuplicate(std::vector<double>& F);                           // Function to check for a non-dominated solution with the same objectives
        int dominion_status(std::vector<double>& F1, std::vector<double>& F2); // 1: F1 weakly dominates F2, 2: F2 dominates F1, 0: neither
        Solution SelectFrom(void);                                          // Function to select a non-dominated solution using the selection metric
        void Calc_Z(void);                                                  // Function to compute the selection metric of the non-dominated solutions
        void SortArchive(std::vector<int>& order, int whichObj);            // Function to sort the non-dominated solutions by an objective
        double neigh_value_continuous(double s, double s_min, double s_max, double r); // Function to perturb a single decision variable

        // Hypervolume functions (Beume and Rudolph, 2006)
        double HV(int data_n, int dim_n, double* ref, double** points);
        void SortPoints(double** X, int size, int which);
        bool covers(double* cub, double* regLow);
        bool partCovers(double* cub, double* regUp);
        int containsBoundary(double* cub, double* regLow, int split);
        double getMeasure(double* regLow, double* regUp);
        int isPile(double* cub, double* regLow, double* regUp);
        double getMedian(double* bounds, int size);
        double computeTrellis(double* regLow, double* regUp, double* trellis);
        void stream(double* regionLow, double* regionUp, double** points, int npoints, int split, double cover);
        int bool_vec_to_ulong(bool* pB, int size);
        void ulong_to_bool_vec(int val, bool* pB, int size);

        // Output functions
        void WriteHeader(FILE* pFile, const char* pAfter);                  // Function to write the objective and parameter column names
        void WriteSolution(FILE* pFile, Solution& solution);                // Function to write the objectives and parameters of a solution
        void WriteArchiveRecord(int iter, double remaining);                // Function to write the non-dominated archive to the run record
        void WriteArchives(FILE* pFile);                                    // Function to write the final non-dominated and dominated archives
        const char* GetSelectionMetricString(void);                         // Function to get the description of the selection metric
};

extern "C" {
void PADDS_Program(int argC, StringType argV[]);
}

#endif /* PADDS_ALGORITHM_H */
