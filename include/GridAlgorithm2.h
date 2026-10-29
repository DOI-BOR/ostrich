/**************************************************************************************************************************************************************

This algorithm will "grid" the objective function surface by evaluating an exhaustive (or nearly exhaustive) set of parameter configurations. Parameter
sets are generated deterministically across a user-specified (or default) number of grid points per parameter, and are dispatched for evaluation in
batches ("mini-grids") that are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

#ifndef GRID_ALGORITHM_H
#define GRID_ALGORITHM_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

// Include C++ classes
#include <algorithm>
#include <vector>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "StatUtility.h"
#include "ParameterGroup.h"
#include "ReadUtility.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"



class GridAlgorithm : public Algorithm {
   public:
        // Define constructors and destructors
        GridAlgorithm();
        ~GridAlgorithm(void){ DBG_PRINT("GridAlgorithm::DTOR"); Destroy(); }
        void Destroy(void);

        // Warm restart is not yet supported for the Grid algorithm
        void WarmStart(void);

        // Solution functions
        void Optimize(void);                                                // Function to perform the exhaustive grid search

        // Output functions
        void WriteStartingMetrics(void);                                     // Function to write the starting Grid metrics
        void WriteEndingMetrics(void);                                       // Function to write the ending Grid metrics

   private:
        // Configuration variables
        std::vector<int> m_GridDims;                                        // Number of grid points along each parameter dimension
        long long m_MiniSize = 0;                                           // Number of grid points evaluated per batch/iteration (0 == use default)

        // Grid geometry, derived from the parameter bounds and m_GridDims
        std::vector<double> m_Lwr;                                         // Lower bound (transformed) for each parameter
        std::vector<double> m_StepSize;                                    // Grid step size (transformed) for each parameter
        std::vector<long long> m_Rval;                                     // Rollover counts used to map a flat grid index to per-parameter indices
        long long m_GridSize = 0;                                          // Total number of points in the exhaustive grid

        // Solution tracking variables
        int m_NumIters = 0;                                                 // Total number of batches/iterations required to cover the grid
        int m_CurIter = 0;                                                  // Counter for the current batch/iteration number

        // Solution functions
        double GetGridVal(long long gridIndex, int paramIndex);            // Computes the value of a parameter at a given flat grid index

};

extern "C" {
void GRID_Program(int argC, StringType argV[]);
}

#endif
