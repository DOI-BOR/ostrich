/**************************************************************************************************************************************************************

Powell's method is a zero-order (derivative-free) optimization algorithm, which utilizes the concept of conjugate directions to determine the optimal search
direction. Each iteration performs a line search along each of the current search directions, followed by a line search along the constructed conjugate
direction. Line searches are performed by the OptSearchClass2 class, with model evaluations solved using the shared, asynchronous MPI work-farming
functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

#ifndef POWELL_ALGORITHM_H
#define POWELL_ALGORITHM_H

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
#include "OptSearchClass2.h"
#include "ParameterGroup.h"
#include "ReadUtility.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class PowellAlgorithm : public Algorithm {
    public:
        // Define constructors and destructors
        PowellAlgorithm();
        ~PowellAlgorithm(void) { DBG_PRINT("PowellAlgorithm::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using Powell's method

        // Output functions
        void WriteStartingMetrics(void);                                    // Function to write the starting Powell metrics
        void WriteEndingMetrics(void);                                      // Function to write the ending Powell metrics

    private:
        // Configuration variables
        int m_MaxIter = 20;                                                 // Maximum number of iterations
        double m_ConvVal = 1E-6;                                            // Convergence value for the change in objective between iterations
        const int m_MaxStallCount = 3;                                      // Number of converged iterations required to stop the analysis

        // Solution tracking variables
        int m_NumDirs = 0;                                                  // Number of search directions
        std::vector<std::vector<double>> m_SearchDirs;                      // Search directions
        OptSearchClass2* m_pSearch = NULL;                                  // 1-dimensional search
        int m_CurIter = 0;                                                  // Current iteration

        // Metrics
        int m_AlgCount = 0;                                                 // Number of evaluations made by the algorithm (outside of the line search)
        int m_NumRestarts = 0;                                              // Number of times the search directions were reset
        int m_NumUprViols = 0;                                              // Counter for upper bound violations
        int m_NumLwrViols = 0;                                              // Counter for lower bound violations

        // Solution functions
        std::vector<double> Evaluate(std::vector<std::vector<double>> samples); // Function to solve a set of parameter sets and track the best solution
        bool LineMinimize(std::vector<double>& X, double& Fcur, std::vector<double>& Scur); // Function to minimize along a direction and step the solution
        void ResetSearchDirections(void);                                   // Function to reset the search directions to the unit vectors
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group
};

extern "C" {
void PWL_Program(int argC, StringType argV[]);
}

#endif /* POWELL_ALGORITHM_H */
