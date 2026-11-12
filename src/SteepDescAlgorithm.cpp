/******************************************************************************
File      : SteepDescAlgorithm.cpp
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of the Steepest Descent optimization algorithm.

The Steepest Descent algorithm is a first-order optimization algorithm which
utilizes the negative of the gradient to determine the search direction. This
implementation has been reworked to utilize the common Algorithm base class
and the new parallel (MPI) evaluation structure that is also used by
GeneticAlgorithm and LevenbergAlgorithm.

Each iteration proceeds in two parallel batch phases:
  1. Gradient phase: the gradient is estimated via finite differences, with
     all of the perturbed parameter sets solved in a single parallel batch
     (reusing the same CalculateJacobianParameters()/CalculateJacobian2()
     helpers used by LevenbergAlgorithm).
  2. Line search phase: rather than the legacy serial Golden-Section/Brent
     search (which evaluates one trial point at a time and therefore cannot
     make use of more than one MPI worker), a batch of candidate step sizes
     is evaluated in parallel each round. If the best candidate lies at the
     edge of the tried range, the batch is regenerated further in that
     direction (a pattern directly analogous to the parallel lambda-trial
     search used by LevenbergAlgorithm); otherwise the interior optimum is
     accepted and the outer iteration completes.

Version History
08-26-03    lsm   created 
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
XX-XX-XX    lsm   Reworked to use the common Algorithm class and the new
                  parallel (MPI) evaluation structure.
******************************************************************************/

// Include the class header
#include "SteepDescAlgorithm.h"

/******************************************************************************
CTOR()

Initializes parameters to reasonable defaults, then overrides them with the
contents of the BeginSteepDescAlg/EndSteepDescAlg configuration block, if
present.
******************************************************************************/
SteepDescAlgorithm::SteepDescAlgorithm(void) {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    //read in algorithm parameters
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        if (CheckToken(pFile, "BeginSteepDescAlg", inputFilename) == true) {
            FindToken(pFile, "EndSteepDescAlg", inputFilename);
            rewind(pFile);

            FindToken(pFile, "BeginSteepDescAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);
            while (strstr(line, "EndSteepDescAlg") == NULL) {
                if (strstr(line, "ConvergenceVal") != NULL) {
                    // Convergence value: searches halt when the relative reduction over
                    // three successive iterations is less than this value
                    sscanf(line, "%s %lf", tmp, &m_ConvVal);

                } else if (strstr(line, "MaxIterations") != NULL) {
                    // Maximum number of iterations in the overall method
                    sscanf(line, "%s %d", tmp, &m_MaxIter);

                } else if (strstr(line, "NumSearchPoints") != NULL) {
                    // Number of candidate step sizes to evaluate in parallel per line-search round.
                    // If unspecified (or <= 0), the number of available MPI workers is used instead.
                    sscanf(line, "%s %d", tmp, &m_NumSearchPoints);

                } else {
                    sprintf(tmp, "Unknown token: %s", line);
                    LogError(ERR_FILE_IO, tmp);
                }

                line = GetNxtDataLine(pFile, inputFilename);
            }
        } else {
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    } else {
        LogError(ERR_FILE_IO, "Using default algorithm setup.");
    }

    IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()

Frees up memory used by the algorithm.
******************************************************************************/
void SteepDescAlgorithm::Destroy(void) {
    IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void SteepDescAlgorithm::WarmStart(void) {
    // todo: warm start has not yet been implemented for algorithms that use the new
    // parallel batch line-search structure (see LevenbergAlgorithm::WarmStart, which
    // has the same limitation on this branch).
}/* end WarmStart() */

/******************************************************************************
CreateInitialStepMultipliers()

Builds the initial batch of candidate step-size multipliers to try along the
current search direction. Multipliers are a symmetric, geometric (base-2)
spread around the current step scale (m_CurStepScale), e.g. for 5 points:
   { scale/4, scale/2, scale, scale*2, scale*4 }

The batch size defaults to the number of available MPI workers (so that a
full round of trials completes in a single parallel pass), but can be
overridden via the NumSearchPoints configuration token.
******************************************************************************/
std::vector<double> SteepDescAlgorithm::CreateInitialStepMultipliers(void) {

    int numPoints = m_NumSearchPoints;
    if (numPoints <= 0) {
        MPI_Comm_size(MPI_COMM_WORLD, &numPoints);
        if (!m_bSolveOnPrimary) { numPoints--; }
    }
    if (numPoints < 3) { numPoints = 3; } //need at least 3 points to distinguish an edge from an interior minimum

    std::vector<double> multipliers;
    int half = numPoints / 2;
    for (int entry = 0; entry < numPoints; entry++) {
        int exponent = entry - half;
        multipliers.push_back(m_CurStepScale * pow(2.0, (double)exponent));
    }

    return multipliers;
}/* end CreateInitialStepMultipliers() */

/******************************************************************************
CreateAdditionalStepMultipliers()

Builds a follow-up batch of candidate step-size multipliers that extends
beyond the edge of the previous batch, in the given direction (+1 to grow the
step further, -1 to shrink it). This mirrors the parallel expansion pattern
used by LevenbergAlgorithm::CreateAdditionalLambdas().
******************************************************************************/
std::vector<double> SteepDescAlgorithm::CreateAdditionalStepMultipliers(double edgeMultiplier, int direction) {

    int numPoints = m_NumSearchPoints;
    if (numPoints <= 0) {
        MPI_Comm_size(MPI_COMM_WORLD, &numPoints);
        if (!m_bSolveOnPrimary) { numPoints--; }
    }
    if (numPoints < 1) { numPoints = 1; }

    std::vector<double> multipliers;
    for (int entry = 1; entry <= numPoints; entry++) {
        multipliers.push_back(edgeMultiplier * pow(2.0, direction * (double)entry));
    }

    return multipliers;
}/* end CreateAdditionalStepMultipliers() */

/******************************************************************************
CreateLineSearchSamples()

Given a base point and a (normalized) search direction, constructs one
candidate parameter set per step-size multiplier: base + multiplier*dir.
Candidates that would fall outside of the parameter bounds are pulled back to
the midpoint between the base value and the violated bound (matching the
boundary-handling behavior of the legacy implementation), and the
corresponding violation counter is incremented.
******************************************************************************/
std::vector<std::vector<double>> SteepDescAlgorithm::CreateLineSearchSamples(std::vector<double> const & basePoint, std::vector<double> const & searchDir,
                                                                              std::vector<double> const & lowerValues, std::vector<double> const & upperValues,
                                                                              std::vector<double> const & multipliers) {
    std::vector<std::vector<double>> samples;

    for (int entryMultiplier = 0; entryMultiplier < multipliers.size(); entryMultiplier++) {
        std::vector<double> candidate = basePoint;

        for (int entryParameter = 0; entryParameter < basePoint.size(); entryParameter++) {
            double tst = basePoint[entryParameter] + (multipliers[entryMultiplier] * searchDir[entryParameter]);

            if (tst > upperValues[entryParameter]) {
                tst = (basePoint[entryParameter] + upperValues[entryParameter]) / 2.00;
                m_NumUprViols++;
            }
            if (tst < lowerValues[entryParameter]) {
                tst = (basePoint[entryParameter] + lowerValues[entryParameter]) / 2.00;
                m_NumLwrViols++;
            }

            candidate[entryParameter] = tst;
        }

        samples.push_back(candidate);
    }

    return samples;
}/* end CreateLineSearchSamples() */

/******************************************************************************
Optimize()

Optimize the objective function using the Steepest Descent algorithm.
******************************************************************************/
void SteepDescAlgorithm::Optimize(void) {

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Perform initial warm start recovery, if needed
    std::vector<int> solveCounts;

    if (m_bWarmStart) {
        std::cout << "Warm start has not been configured for the Steepest Descent Algorithm. Exiting the analysis..." << std::endl;
        throw std::invalid_argument("Warm start has not been configured for the Steepest Descent Algorithm");

    } else {
        // Fill the solve counter with zeros
        for (int entryWorker = 0; entryWorker < numberOfMpiProcesses; entryWorker++) {
            solveCounts.push_back(0);
        }
    }

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    // Write setup and starting metrics
    WriteSetup2("SteepestDescent", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), m_pParamGroup->GetNumParams(), m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();

    // Get the initial parameter values and bounds out of the parameter group
    std::vector<double> lowerValues, upperValues;
    std::vector<double> currentValues;

    for (int entry = 0; entry < m_pParamGroup->GetNumParams(); entry++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entry);

        currentValues.push_back(temp->GetEstimatedValueTransformed());
        upperValues.push_back(temp->GetUpperBoundTransformed());
        lowerValues.push_back(temp->GetLowerBoundTransformed());
    }

    // Evaluate the initial (starting) point
    std::vector<std::vector<double>> initialSample = { currentValues };
    std::vector<double> initialObjective(1, INFINITY);
    ManageSingleObjectiveIterations(initialSample, 0, m_pParamGroup->GetNumParams(), initialObjective);

    m_BestObjective = initialObjective[0];
    m_BestAlternative = currentValues;

    // Write the starting (iteration 0) record
    WriteRecord2(0, m_BestObjective, 0.0, m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Iteration bookkeeping
    m_Iteration = 0;
    m_dObjFuncHistory.clear();
    bool continueIterations = true;

    while (continueIterations) {

        m_Iteration++;
        m_CurIter = m_Iteration;

        // ---------------------------------------------------------------
        // Gradient phase: estimate the gradient using a single parallel
        // batch of finite-difference evaluations.
        // ---------------------------------------------------------------
        std::vector<std::vector<double>> jacobianSamples;
        std::vector<bool> lockedParameters;
        CalculateJacobianParameters(m_BestAlternative, lowerValues, upperValues, m_StepSize, m_StepToleranceMaximum, m_DerivativeType, jacobianSamples, lockedParameters);

        std::vector<double> jacobianObjectives(jacobianSamples.size(), INFINITY);
        ManageSingleObjectiveIterations(jacobianSamples, 0, m_pParamGroup->GetNumParams(), jacobianObjectives);

        std::vector<std::vector<double>> jacobian;
        CalculateJacobian2(m_BestAlternative, m_BestObjective, jacobianObjectives, jacobianSamples, m_DerivativeType, lockedParameters, jacobian);

        std::vector<double> gradient = jacobian.back();

        // Determine the largest-magnitude gradient component
        double maxAbsGrad = 0.00;
        for (int entry = 0; entry < gradient.size(); entry++) {
            double mag = fabs(gradient[entry]);
            if (mag > maxAbsGrad) { maxAbsGrad = mag; }
        }

        double bestLineObjective = m_BestObjective;
        std::vector<double> bestLineAlternative = m_BestAlternative;

        if (maxAbsGrad > NEARLY_ZERO) {
            // -------------------------------------------------------------
            // Build the normalized descent direction (negative gradient,
            // scaled so that the largest-magnitude component is 1.00)
            // -------------------------------------------------------------
            std::vector<double> searchDir(gradient.size());
            for (int entry = 0; entry < gradient.size(); entry++) {
                searchDir[entry] = (-1.00 * gradient[entry]) / maxAbsGrad;
            }

            // Initialize the line-search step scale, if not already set
            if (m_CurStepScale <= 0.00) {
                double sumRange = 0.00;
                int rangeCount = 0;
                for (int entry = 0; entry < lockedParameters.size(); entry++) {
                    if (!lockedParameters[entry]) {
                        sumRange += (upperValues[entry] - lowerValues[entry]);
                        rangeCount++;
                    }
                }
                m_CurStepScale = (rangeCount > 0) ? (0.05 * (sumRange / rangeCount)) : 0.05;
            }

            // -------------------------------------------------------------
            // Line search phase: evaluate a batch of candidate step sizes
            // in parallel, extending the batch toward the edge that
            // improves upon the current best, until an interior optimum is
            // found (or no improving direction can be found).
            // -------------------------------------------------------------
            double bestMultiplierUsed = m_CurStepScale;
            std::vector<double> multipliers = CreateInitialStepMultipliers();

            const int maxLineSearchRounds = 8;
            int lineSearchRound = 0;
            bool lineSearchDone = false;

            while (!lineSearchDone) {
                std::vector<std::vector<double>> lineSamples = CreateLineSearchSamples(m_BestAlternative, searchDir, lowerValues, upperValues, multipliers);

                std::vector<double> lineObjectives(lineSamples.size(), INFINITY);
                ManageSingleObjectiveIterations(lineSamples, 0, m_pParamGroup->GetNumParams(), lineObjectives);

                double roundBestObjective;
                int roundBestIndex;
                GetBestSingleObjective(lineObjectives, roundBestObjective, roundBestIndex);

                if (roundBestObjective < bestLineObjective) {
                    bestLineObjective = roundBestObjective;
                    bestLineAlternative = lineSamples[roundBestIndex];
                    bestMultiplierUsed = multipliers[roundBestIndex];
                }

                lineSearchRound++;

                if (multipliers.size() > 1 && lineSearchRound < maxLineSearchRounds) {
                    int minIdx = std::min_element(multipliers.begin(), multipliers.end()) - multipliers.begin();
                    int maxIdx = std::max_element(multipliers.begin(), multipliers.end()) - multipliers.begin();

                    if (roundBestObjective >= m_BestObjective) {
                        // No improvement found this round.
                        if (roundBestIndex == minIdx) {
                            // The least-bad trial is at the smallest step tried; try smaller steps once more.
                            multipliers = CreateAdditionalStepMultipliers(multipliers[minIdx], -1);
                        } else {
                            // No improving direction found at this scale; stop searching.
                            lineSearchDone = true;
                        }
                    } else if (roundBestIndex == maxIdx) {
                        // Still improving at the largest step tried; extend further.
                        multipliers = CreateAdditionalStepMultipliers(multipliers[maxIdx], 1);
                    } else if (roundBestIndex == minIdx) {
                        // Improvement found only at the smallest step tried; try smaller steps.
                        multipliers = CreateAdditionalStepMultipliers(multipliers[minIdx], -1);
                    } else {
                        // Improvement found at an interior point; accept.
                        lineSearchDone = true;
                    }
                } else {
                    lineSearchDone = true;
                }
            }

            // Adapt the step scale for the next iteration's line search
            if (bestLineObjective < m_BestObjective) {
                m_CurStepScale = bestMultiplierUsed;
            } else {
                m_CurStepScale = m_CurStepScale / 4.00;
            }
        }/* end if(gradient is non-zero) */

        // Compute the relative reduction for this iteration
        double dObjFunc = fabs((m_BestObjective - bestLineObjective) / (m_BestObjective + NEARLY_ZERO));

        // Update the best solution, if improved
        if (bestLineObjective < m_BestObjective) {
            m_BestObjective = bestLineObjective;
            m_BestAlternative = bestLineAlternative;
        }

        // Set the best values into the parameter group (for reporting/output)
        for (int entryParam = 0; entryParam < m_pParamGroup->GetNumParams(); entryParam++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParam);
            temp->SetEstimatedValueTransformed(m_BestAlternative[entryParam]);
        }

        // Write the iteration result to the log file
        WriteRecord2(m_Iteration, m_BestObjective, dObjFunc, m_pObsGroup, m_pObjFunc, m_pParamGroup);

        // Track relative-reduction history: convergence requires the reduction over the
        // last three successive iterations to fall below the configured ConvergenceVal.
        m_dObjFuncHistory.push_back(dObjFunc);
        if (m_dObjFuncHistory.size() > 3) {
            m_dObjFuncHistory.erase(m_dObjFuncHistory.begin());
        }

        bool converged = false;
        if (m_dObjFuncHistory.size() == 3) {
            double maxRecentChange = *std::max_element(m_dObjFuncHistory.begin(), m_dObjFuncHistory.end());
            if (maxRecentChange < m_ConvVal) { converged = true; }
        }

        if (converged || (m_Iteration >= m_MaxIter) || (maxAbsGrad <= NEARLY_ZERO)) {
            continueIterations = false;
        }
    }/* end while(continueIterations) */

    // Terminate the secondary workers and let them clean up.
    TerminateWorkers();

    //write algorithm metrics
    WriteEndingMetrics();
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using the Steepest Descent algorithm.

Note: post-hoc regression statistics (variance/covariance, via StatsClass)
have not yet been ported to work with the new Algorithm base class - the same
limitation currently applies to GeneticAlgorithm::Calibrate() and
LevenbergAlgorithm::Calibrate() on this branch. For now, Calibrate() simply
performs the optimization.
******************************************************************************/
void SteepDescAlgorithm::Calibrate(void) {
    Optimize();
}/* end Calibrate() */

/******************************************************************************
WriteStartingMetrics()

Write out algorithm setup information prior to running the analysis.
******************************************************************************/
void SteepDescAlgorithm::WriteStartingMetrics(void) {

    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    fprintf(pFile, "Algorithm         : Steepest-Descent\n");
    fprintf(pFile, "Max Iterations    : %d\n", m_MaxIter);
    fprintf(pFile, "Convergence Val   : %E\n", m_ConvVal);
    fprintf(pFile, "Search Points     : %s\n", (m_NumSearchPoints > 0) ? std::to_string(m_NumSearchPoints).c_str() : "auto (# of MPI workers)");

    // Write the parameter banner to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Convergence\n");

    fclose(pFile);
}/* end WriteStartingMetrics() */

/******************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
******************************************************************************/
void SteepDescAlgorithm::WriteEndingMetrics(void) {

    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm        : Steepest-Descent\n");
    fprintf(pFile, "Max Iterations   : %d\n", m_MaxIter);
    fprintf(pFile, "Convergence Val  : %E\n", m_ConvVal);
    fprintf(pFile, "Iterations       : %d\n", m_CurIter);
    fprintf(pFile, "Total Evals      : %d\n", m_NumSolves);
    fprintf(pFile, "Upper Violations : %d\n", m_NumUprViols);
    fprintf(pFile, "Lower Violations : %d\n", m_NumLwrViols);

    fclose(pFile);
}/* end WriteEndingMetrics() */

/******************************************************************************
STPDSC_Program()

Create a model and solve using the steepest descent algorithm.
******************************************************************************/
void STPDSC_Program(int argc, StringType argv[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine the behavior based on rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine.
        SteepDescAlgorithm* SteepDesc = new SteepDescAlgorithm();

        // Begin the primary worker operations. This transfers information to the secondary
        // workers and manages the solution process.
        SteepDesc->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine.
        ModelWorker worker = ModelWorker(true);

        // Start the work on the processor.
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions.
    }

} /* end STPDSC_Program() */
