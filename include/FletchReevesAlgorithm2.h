/**************************************************************************************************************************************************************

The Fletcher-Reeves algorithm is a first-order optimization algorithm, which utilizes the concept of conjugate directions in conjunction with the steepest-
descent information (negative of the gradient). Each iteration performs a line search along the current conjugate direction, computes the gradient at the
new location, and updates the conjugate direction using the Polak-Ribiere formula (reverting to steepest descent if the direction is not a descent
direction). The gradient is computed by the OptMathClass2 class and the line searches are performed by the OptSearchClass2 class, with model evaluations
solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

#ifndef FLETCH_REEVES_ALGORITHM_H
#define FLETCH_REEVES_ALGORITHM_H

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
#include "OptMathClass2.h"
#include "OptSearchClass2.h"
#include "ParameterGroup.h"
#include "ReadUtility.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class FletchReevesAlgorithm : public Algorithm {
    public:
        // Define constructors and destructors
        FletchReevesAlgorithm();
        ~FletchReevesAlgorithm(void) { DBG_PRINT("FletchReevesAlgorithm::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using Fletcher-Reeves

        // Output functions
        void WriteStartingMetrics(void);                                    // Function to write the starting Fletcher-Reeves metrics
        void WriteEndingMetrics(void);                                      // Function to write the ending Fletcher-Reeves metrics

    private:
        // Configuration variables
        int m_MaxIter = 20;                                                 // Maximum number of iterations
        double m_ConvVal = 1E-6;                                            // Convergence value for the change in objective between iterations
        int m_MaxCount = 3;                                                 // Number of converged iterations required to stop the analysis

        // Solution tracking variables
        int m_NumParams = 0;                                                // Number of parameters
        OptMathClass2* m_pMath = NULL;                                      // Gradient calculations
        OptSearchClass2* m_pSearch = NULL;                                  // 1-dimensional search
        int m_CurIter = 0;                                                  // Current iteration

        // Metrics
        int m_AlgCount = 0;                                                 // Number of evaluations made by the algorithm (outside of the gradient and line search)
        int m_NumRestarts = 0;                                              // Number of times the search reverted to steepest descent
        int m_NumUprViols = 0;                                              // Counter for upper bound violations
        int m_NumLwrViols = 0;                                              // Counter for lower bound violations

        // Solution functions
        std::vector<double> Evaluate(std::vector<std::vector<double>> samples); // Function to solve a set of parameter sets and track the best solution
        void SetParameterGroup(std::vector<double> parameters);             // Function to set values into the parameter group
};

extern "C" {
void FLRV_Program(int argC, StringType argV[]);
}

#endif /* FLETCH_REEVES_ALGORITHM_H */
