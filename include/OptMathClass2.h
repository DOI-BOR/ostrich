/**************************************************************************************************************************************************************

Computes the gradient of the objective function using forward finite differences, for use by gradient-based optimization algorithms. Unlike the pre-refactor
OptMathClass, this class does not interact with a model directly. Instead, parameter sets are evaluated through a user-supplied evaluator, which allows the
perturbations to be solved concurrently using the asynchronous MPI work-farming functionality provided by the Algorithm base class.

The finite difference increments are configured using the optional BeginMathAndStats group (DiffIncType, DiffRelIncrement, DiffIncrement and
DiffMinIncrement).

**************************************************************************************************************************************************************/

#ifndef OPT_MATH_CLASS2_H
#define OPT_MATH_CLASS2_H

// Include C classes
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Include C++ classes
#include <functional>
#include <vector>

// Include custom classes
#include "MyHeaderInc.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class OptMathClass2 {
    public:
        // Function used to evaluate a set of parameter alternatives and return their objectives (in the same order)
        typedef std::function<std::vector<double>(std::vector<std::vector<double>>)> Evaluator;

        // Define constructors and destructors
        OptMathClass2(ParameterGroup* pParamGroup, Evaluator evaluator);
        ~OptMathClass2(void) { DBG_PRINT("OptMathClass2::DTOR"); Destroy(); }
        void Destroy(void);

        // Solution functions
        std::vector<double> CalcGradient(std::vector<double> point, double objective, double& fmin,
                                         std::vector<double>& xmin);  // Function to compute the gradient at a point

        // Output functions
        void WriteMetrics(FILE* pFile);                                // Function to write the finite difference setup and metrics

    private:
        // Configuration variables
        FiniteDiffIncType m_DiffIncType = FD_RANGE_REL;                // Increment type used to compute the derivatives
        std::vector<double> m_DiffInc;                                 // Finite difference increment for each parameter
        double m_MinInc = NEARLY_ZERO;                                 // Smallest allowable (value-relative) increment

        // Solution tracking variables
        ParameterGroup* m_pParamGroup = NULL;                          // Parameter group (used for bounds)
        Evaluator m_Evaluate;                                          // Function used to evaluate parameter sets
        int m_NumParams = 0;                                           // Number of parameters

        // Metrics
        int m_GradCount = 0;                                           // Number of gradient calculations
        int m_DiffCount = 0;                                           // Number of evaluations used to compute derivatives

        // Solution functions
        void InitFromFile(IroncladString pFileName);                   // Function to read the finite difference configuration
        double CalcIncrement(int index, double value, FiniteDiffIncType incrementType); // Function to compute the signed increment for a parameter
        void CalcDerivatives(std::vector<int> indices, std::vector<double> point, double objective, FiniteDiffIncType incrementType,
                             std::vector<double>& gradient, double& fmin, std::vector<double>& xmin); // Function to compute a set of derivatives together
};

#endif /* OPT_MATH_CLASS2_H */
