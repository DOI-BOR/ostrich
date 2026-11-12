/******************************************************************************
File      : SteepDescAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of the Steepest Descent optimization algorithm.

Version History
08-26-03    lsm   created 
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
XX-XX-XX    lsm   Reworked to use the common Algorithm class and the new
                  parallel (MPI) evaluation structure. The 1-D line search
                  is now performed by evaluating a batch of candidate step
                  sizes in parallel each round (similar in spirit to the
                  parallel lambda-trial search used by LevenbergAlgorithm),
                  rather than the old serial Golden-Section/Brent search.
******************************************************************************/
#ifndef STEEP_DESC_ALGORITHM_H
#define STEEP_DESC_ALGORITHM_H

#include "MyHeaderInc.h"

// Include C classes
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

// Include custom classes
#include "Algorithm.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "MathClass2.h"
#include "Exception.h"
#include "Utility.h"
#include "WriteUtility2.h"

/******************************************************************************
class SteepDescAlgorithm

The Steepest Descent algorithm is a first-order optimization algorithm, which 
utilizes the negative of the gradient to determine the optimal search 
direction. Each iteration computes the gradient (via parallel finite-
difference evaluations) and then performs a parallel batch line search along 
the resulting descent direction to determine the step size.
******************************************************************************/
class SteepDescAlgorithm : public Algorithm
{
   public:
      SteepDescAlgorithm(void);
      ~SteepDescAlgorithm(void){ DBG_PRINT("SteepDescAlgorithm::DTOR"); Destroy(); }
      void Destroy(void);

      void WarmStart(void);

      // Solution functions
      void Optimize(void);
      void Calibrate(void);

      // Output functions
      void WriteStartingMetrics(void);
      void WriteEndingMetrics(void);

   private:
      // Helper routines for the parallel batch line search
      std::vector<double> CreateInitialStepMultipliers(void);
      std::vector<double> CreateAdditionalStepMultipliers(double bestMultiplier, int direction);
      std::vector<std::vector<double>> CreateLineSearchSamples(std::vector<double> const & basePoint, std::vector<double> const & searchDir,
                                                                std::vector<double> const & lowerValues, std::vector<double> const & upperValues,
                                                                std::vector<double> const & multipliers);

      //Configuration variables (read from the BeginSteepDescAlg block)
      double m_ConvVal = 1E-6;      //convergence value: searches halt when the relative reduction over
                                     //three successive iterations is less than this value
      int m_MaxIter = 20;           //maximum number of iterations in the overall method
      int m_NumSearchPoints = 0;    //number of candidate step sizes to try in parallel per line-search
                                     //round (0 = size to the number of available MPI workers)

      //finite-difference settings used to compute the gradient (kept as internal
      //defaults, consistent with the legacy OptMathClass behavior)
      double m_StepSize = 0.001;         //relative finite-difference step used to compute the gradient
      double m_StepToleranceMaximum = 0; //parameters whose perturbed range falls below this are locked
      int m_DerivativeType = FIRST_FORWARD;

      //adaptive scale (in normalized search-direction units) used to seed each iteration's
      //parallel batch line search; initialized on first use, then adapted iteration-to-iteration
      double m_CurStepScale = -1.00;

      //Solution tracking variables
      int m_Iteration = 0;
      int m_CurIter = 0;

      //metrics
      int m_AlgCount = 0;
      int m_NumUprViols = 0;
      int m_NumLwrViols = 0;

      //relative-reduction history used for the 3-iteration convergence check
      std::vector<double> m_dObjFuncHistory;
}; /* end class SteepDescAlgorithm */

extern "C" {
void STPDSC_Program(int argc, StringType argv[]);
}

#endif /* STEEP_DESC_ALGORITHM_H */

