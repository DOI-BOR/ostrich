/**************************************************************************************************************************************************************

One-dimensional line search used to compute the optimum step size along a given search direction. The minimum is first bracketed and then refined using
either the Golden Section or Brent method. Unlike the pre-refactor OptSearchClass, this class does not interact with a model directly. Instead, parameter sets
are evaluated through a user-supplied evaluator, which allows the search to be driven by the asynchronous MPI work-farming functionality provided by the
Algorithm base class. Independent evaluations are submitted to the evaluator together so that they can be solved concurrently.

**************************************************************************************************************************************************************/

#ifndef OPT_SEARCH_CLASS2_H
#define OPT_SEARCH_CLASS2_H

// Include C classes
#include <math.h>
#include <stdio.h>
#include <string.h>

// Include C++ classes
#include <functional>
#include <vector>

// Include custom classes
#include "MyHeaderInc.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "WriteUtility2.h"
#include "Exception.h"
#include "Utility.h"


class OptSearchClass2 {
    public:
        // Function used to evaluate a set of parameter alternatives and return their objectives (in the same order)
        typedef std::function<std::vector<double>(std::vector<std::vector<double>>)> Evaluator;

        // Available 1-D search methods
        enum SearchMethod { SEARCH_GOLDEN_SECTION = 0, SEARCH_BRENT = 1 };

        // Define constructors and destructors
        OptSearchClass2(ParameterGroup* pParamGroup, Evaluator evaluator);
        ~OptSearchClass2(void) { DBG_PRINT("OptSearchClass2::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        double CalcStepSize(std::vector<double> startPoint, double startObjective, std::vector<double> direction, double& fmin,
                            std::vector<double>& xmin);          // Function to compute the optimum step size along a direction

        // Output functions
        void WriteMetrics(FILE* pFile);                          // Function to write the line search setup and metrics

    private:
        // Bracketing triplet of step sizes and associated objectives
        struct Bracket {
            double a, b, c;
            double fa, fb, fc;
        };

        // Step size and associated objective
        struct MinimumPoint {
            double x, fx;
        };

        // Configuration variables
        double m_SearchConvVal = 1E-4;                           // Convergence criteria for the 1-D search
        int m_SearchType = SEARCH_GOLDEN_SECTION;                // Type of 1-D search
        const int m_MaxBrentIters = 36;                          // Maximum number of Brent iterations before reverting to Golden Section

        // Solution tracking variables
        ParameterGroup* m_pParamGroup = NULL;                    // Parameter group (used for bounds)
        Evaluator m_Evaluate;                                    // Function used to evaluate parameter sets
        int m_NumParams = 0;                                     // Number of parameters
        std::vector<double> m_StartPoint;                        // Location the search is started from
        double m_StartObjective = INFINITY;                      // Objective at the start location
        std::vector<double> m_Dir;                               // Search direction

        // Metrics
        int m_BoundMinCount = 0;                                 // Evaluations used to bracket the minimum
        int m_GoldSectCount = 0;                                 // Evaluations used by the Golden Section method
        int m_BrentCount = 0;                                    // Evaluations used by the Brent method

        // Solution functions
        void InitFromFile(IroncladString pFileName);             // Function to read the 1-D search configuration
        std::vector<double> CreateAlphaPoint(double alpha);      // Function to create the parameter set for a step size
        std::vector<double> CalcF(std::vector<double> alphas, double& fmin, std::vector<double>& xmin); // Function to evaluate a set of step sizes
        double CalcF(double alpha, double& fmin, std::vector<double>& xmin); // Function to evaluate a single step size
        double LimitStepSize(double alpha);                      // Function to prevent a step from exceeding the parameter bounds
        Bracket BracketMinimum(double a, double b, double& fmin, std::vector<double>& xmin); // Function to bracket the minimum
        MinimumPoint GoldSect(Bracket bracket, double& fmin, std::vector<double>& xmin); // Function to refine the minimum with Golden Section
        MinimumPoint Brent(Bracket bracket, double& fmin, std::vector<double>& xmin); // Function to refine the minimum with Brent's method
};

#endif /* OPT_SEARCH_CLASS2_H */
