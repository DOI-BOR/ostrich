/**************************************************************************************************************************************************************

The Bisection algorithm optimizes one parameter at a time via simple bisection (a coordinate-descent style search). For each parameter, the algorithm 
subdivides the parameter's bounded range into quadrants, evaluates the objective function at the resulting quadrant boundaries, and repeatedly discards 
the quadrant(s) least likely to contain the optimum. Parameters are processed one at a time (holding all other parameters fixed at their current best 
values) and the full cycle through all parameters is repeated for a configurable number of outer iterations, with a random perturbation of parameter 
values in between the first and second outer iterations onward.

**************************************************************************************************************************************************************/

#ifndef BISECTION_ALGORITHM_H
#define BISECTION_ALGORITHM_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Include C++ classes
#include <vector>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class BisectionAlgorithm : public Algorithm {
   public:
        // Define constructors and destructors
        BisectionAlgorithm();
        ~BisectionAlgorithm(void){ DBG_PRINT("BisectionAlgorithm::DTOR"); Destroy(); }
        void Destroy(void);

        // Define sampling functions used by the Bisection algorithm
        void WarmStart(void);                                              // todo: implement - placeholder for warm starting the search from a previous run

        // Solution functions
        void Optimize(void);                                               // Function to optimize a solution using bisection
        void Calibrate(void);                                              // Function to solve least squares fitting using bisection

        // Output functions
        void WriteStartingMetrics(void);                                   // Function to write the starting Bisection metrics
        void WriteEndingMetrics(void);                                     // Function to write the ending Bisection metrics


   private:
        // Configuration variables
        int m_MaxOuter = 50;                                               // Number of outer iterations (one full cycle through all parameters)
        int m_MaxInner = 20;                                                // Number of inner iterations (bisections) applied to each parameter
        double m_ConvVal = 1E-6;                                            // Convergence value (parsed for parity with the legacy configuration but,
                                                                             // consistent with the pre-existing implementation, not currently used to
                                                                             // trigger early termination)

        // Solution tracking variables
        int m_NumParams = 0;                                                // Number of parameters being optimized
        int m_CurIter = 0;                                                  // Current outer iteration
        int m_AlgCount = 0;                                                 // Number of model evaluations performed by the algorithm

        // Solution functions
        std::vector<double> BuildCandidate(std::vector<double> baseValues, int paramIndex, double paramValue); // Builds a full parameter vector, substituting a single parameter's value

}; 

extern "C" {
void BIS_Program(int argC, StringType argV[]);
}

#endif 
