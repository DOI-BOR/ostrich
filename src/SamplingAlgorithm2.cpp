/**************************************************************************************************************************************************************

An implementation of a sampling algorithm, loosely based on Big Bang-Big Crunch (BB-BC). An initial set of Latin Hypercube samples is evaluated, and each
subsequent iteration samples the space around the current best solution using a search radius that shrinks as the analysis progresses. The samples of each
iteration are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "SamplingAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the sampling algorithm. First, all parameters are assigned default values and then the user input file is checked for overriding values. The
sampling schedule (number of iterations and samples per iteration) is computed from the evaluation budget once it is known.
**************************************************************************************************************************************************************/
SamplingAlgorithm::SamplingAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    int numberOfParameters = m_pParamGroup->GetNumParams();
    IroncladString inputFilename = GetInFileName();

    // Read sampling algorithm information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Sampling Algorithm block within the file
        if (CheckToken(pFile, "BeginSamplingAlg", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndSamplingAlg", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginSamplingAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndSamplingAlg") == NULL) {

                if (strstr(line, "MaxEvaluations") != NULL) {
                    // Set the evaluation budget
                    sscanf(line, "%s %d", tmp, &m_MaxEvals);
                    if (m_MaxEvals < 1) { m_MaxEvals = 1; }
                }

                // Get the next line
                line = GetNxtDataLine(pFile, inputFilename);
            }
        }
        else {
            // The sampling algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Split the evaluation budget into iterations and samples per iteration
    m_MaxIter = (int)(sqrt((double)m_MaxEvals));
    m_NumSamples = m_MaxEvals / m_MaxIter;
    m_NumExtra = m_MaxEvals - (m_MaxIter * m_NumSamples);

    // Initialize the per-parameter bounds, standard deviations and forward weights
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
        double lwr = temp->GetLowerBoundTransformed();
        double upr = temp->GetUpperBoundTransformed();

        m_Lwr.push_back(lwr);
        m_Upr.push_back(upr);
        m_SD.push_back(0.25 * (upr - lwr));
        m_Fwd.push_back(0.50);
    }

    // TODO: port support for user-defined evaluations via a BeginInitParams/EndInitParams block, similar to what the
    // pre-refactor SamplingAlgorithm supported. This has not yet been ported to the new Algorithm-based structure
    // (GeneticAlgorithm and ParticleSwarm have the same limitation currently).

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the sampling algorithm and it's member variables
**************************************************************************************************************************************************************/
void SamplingAlgorithm::Destroy(void) {

    m_Lwr.clear();
    m_Upr.clear();
    m_SD.clear();
    m_Fwd.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
CreateInitialSample()

Creates the initial set of samples using Latin Hypercube sampling across the parameter bounds.
**************************************************************************************************************************************************************/
std::vector<std::vector<double>> SamplingAlgorithm::CreateInitialSample(int sampleSize) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Use Latin Hypercube to initialize the space
    LatinHypercube* pLHS = new LatinHypercube(sampleSize, numberOfParameters, sampleSize);
    MEM_CHECK(pLHS);

    // Fill the array from the LHS object
    pLHS->CreateUniformSample(m_Lwr.data(), m_Upr.data());
    std::vector<std::vector<double>> samples = pLHS->GetSampleMatrix();

    // Cleanup the memory
    delete pLHS;

    // Return to the calling function
    return samples;
}

/**************************************************************************************************************************************************************
CreateSample()

Creates a set of random samples around the current best solution. Each parameter is displaced in a random direction (biased by the forward weights) by a
random fraction of the search radius, which shrinks linearly over the course of the analysis and is limited so that samples remain within the bounds.

The pre-refactor algorithm generated and solved samples one at a time, reversing the forward weights after every sample that did not improve the best
solution. All samples of an iteration are now generated before any are solved, so the weights are reversed between consecutive samples here to reproduce
the sequential behavior for an iteration without improvement.
**************************************************************************************************************************************************************/
std::vector<std::vector<double>> SamplingAlgorithm::CreateSample(int sampleSize) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    std::vector<std::vector<double>> samples;

    for (int entrySample = 0; entrySample < sampleSize; entrySample++) {
        std::vector<double> temp;

        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            double best = m_BestAlternative[entryParameter];
            double fwd = m_Fwd[entryParameter];
            if ((entrySample % 2) == 1) { fwd = 1.00 - fwd; }

            // Compute the direction of displacement
            double sign, maxAlpha;
            double r = (double)MyRand() / (double)MY_RAND_MAX;
            if (r < fwd) {
                sign = +1.00;
                maxAlpha = (m_Upr[entryParameter] - best) / m_SD[entryParameter]; // max positive displacement
            }
            else {
                sign = -1.00;
                maxAlpha = (best - m_Lwr[entryParameter]) / m_SD[entryParameter]; // max negative displacement
            }

            // Compute the magnitude of displacement
            double alpha = m_Radius * ((m_MaxIter - (m_CurIter - 1)) / (double)m_MaxIter);
            if (alpha > maxAlpha) { alpha = maxAlpha; }
            r = (double)MyRand() / (double)MY_RAND_MAX;
            alpha *= r;

            temp.push_back(best + (sign * alpha * m_SD[entryParameter]));
        }

        samples.push_back(temp);
    }

    // Return to the calling function
    return samples;
}

/**************************************************************************************************************************************************************
UpdateBest()

Processes a set of solved samples in order. A sample that improves on the best solution becomes the new best, the standard deviations are rescaled from it,
and the forward weights are pointed in the direction that produced the improvement. A sample that does not improve on the best solution reverses the
forward weights.
**************************************************************************************************************************************************************/
void SamplingAlgorithm::UpdateBest(std::vector<std::vector<double>> samples, std::vector<double> objectives) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    for (int entrySample = 0; entrySample < (int)samples.size(); entrySample++) {
        if (m_BestAlternative.empty() || objectives[entrySample] < m_BestObjective) {
            // The sample is the new best solution
            std::vector<double> previousBest = m_BestAlternative;
            m_BestObjective = objectives[entrySample];
            m_BestAlternative = samples[entrySample];

            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                m_SD[entryParameter] = m_BestAlternative[entryParameter] * ((m_MaxIter - m_CurIter) / (double)m_MaxIter);

                // Point the forward weights in the direction of improvement (unchanged for the first solution)
                if (!previousBest.empty()) {
                    if (m_BestAlternative[entryParameter] < previousBest[entryParameter]) { m_Fwd[entryParameter] = 0.00; }
                    else { m_Fwd[entryParameter] = 1.00; }
                }
            }
        }
        else {
            // Reverse direction
            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                m_Fwd[entryParameter] = 1.00 - m_Fwd[entryParameter];
            }
        }
    }
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using the sampling algorithm.
**************************************************************************************************************************************************************/
void SamplingAlgorithm::Optimize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Perform initial warm start recovery, if needed
    std::vector<int> solveCounts;
    std::vector<std::vector<double>> iterationValues;
    std::vector<std::vector<std::vector<double>>> workerValues;

    if (m_bWarmStart) {
        // Read model file for previous iteration count
        iterationValues = ReadOutputFile();

        // Read worker files for total number of solves
        workerValues = ReadWorkerFiles(numberOfMpiProcesses, m_bSolveOnPrimary);

        // Read model file for previous best
        FindPreviousBest(workerValues);

        // Adjust the primary counts
        AdjustPrimaryCounts(workerValues);

        // Adjust the worker counts
        solveCounts = AdjustSecondaryCounts(workerValues);

        // Update cache to prevent multiple solutions
        UpdateCache(workerValues);
    }
    else {
        // Fill the solve counter with zeros
        solveCounts = std::vector<int>(numberOfMpiProcesses, 0);
    }

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    if (m_bWarmStart && iterationValues.size() > 0 && !m_BestAlternative.empty()) {
        // Resume at the iteration following the last one recorded in the log file
        m_CurIter = (int)iterationValues.size();
        m_AlgCount = m_NumSolves;

        // Rescale the standard deviations from the recovered best solution. The forward weights cannot be recovered and retain their defaults.
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            m_SD[entryParameter] = m_BestAlternative[entryParameter] * ((m_MaxIter - m_CurIter) / (double)m_MaxIter);
        }
    }
    else {
        if (!m_bWarmStart) {
            // Start a clean analysis
            // Write the information to the primary log file
            WriteSetup2("SamplingAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
            WriteStartingMetrics();
        }

        // Generate and evaluate the initial LHS samples
        m_CurIter = 0;
        std::vector<std::vector<double>> samples = CreateInitialSample(m_NumSamples + m_NumExtra);
        std::vector<double> objectives(samples.size(), INFINITY);
        ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);
        m_AlgCount += (int)samples.size();

        UpdateBest(samples, objectives);

        // Set the values into the pointer groups
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
            temp->SetEstimatedValueTransformed(m_BestAlternative[entryParameter]);
        }

        // Write the initial configuration to the log file
        WriteRecord2(m_CurIter, m_BestObjective, 100.00, m_pObsGroup, m_pObjFunc, m_pParamGroup);

        m_CurIter = 1;
    }

    // Main optimization loop. Each iteration samples the space around the current best solution.
    while (m_CurIter < m_MaxIter) {

        double previousBest = m_BestObjective;

        // Generate and evaluate the samples for this iteration
        std::vector<std::vector<double>> samples = CreateSample(m_NumSamples);
        std::vector<double> objectives(samples.size(), INFINITY);
        ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);
        m_AlgCount += (int)samples.size();

        UpdateBest(samples, objectives);

        // Compute the relative change (percent) in the best objective
        double relativeChange = fabs(100.00 * ((previousBest - m_BestObjective) / previousBest));

        // Set the values into the pointer groups
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
            temp->SetEstimatedValueTransformed(m_BestAlternative[entryParameter]);
        }

        // Write the iteration result to the log file
        WriteRecord2(m_CurIter, m_BestObjective, relativeChange, m_pObsGroup, m_pObjFunc, m_pParamGroup);

        m_CurIter++;
    }

    // Terminate the secondary workers and let them clean up
    TerminateWorkers();

    // Write algorithm metrics
    WriteEndingMetrics();

} /* end Optimize() */

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for the sampling algorithm.
**************************************************************************************************************************************************************/
void SamplingAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Evaluations         : %d\n", m_MaxEvals);
    fprintf(pFile, "Max Iterations          : %d\n", m_MaxIter);
    fprintf(pFile, "Samples per Iteration   : %d\n", m_NumSamples);
    fprintf(pFile, "Extra Initial Samples   : %d\n", m_NumExtra);

    // Write the parameter header to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "RelativeChange\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void SamplingAlgorithm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Sampling Method\n");
    fprintf(pFile, "Max Evaluations         : %d\n", m_MaxEvals);
    fprintf(pFile, "Iterations              : %d\n", m_CurIter);
    fprintf(pFile, "Samples per Iteration   : %d\n", m_NumSamples);
    fprintf(pFile, "Extra Samples           : %d\n", m_NumExtra);
    fprintf(pFile, "Algorithm Evals         : %d\n", m_AlgCount);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);

    fprintf(pFile, "\nParameter Standard Deviations (final estimate)\n");
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->Write(pFile, WRITE_BNR);
        fprintf(pFile, " : %E\n", m_SD[entryParameter]);
    }

    fprintf(pFile, "\nParameter Forward Weights (at final iteration)\n");
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->Write(pFile, WRITE_BNR);
        fprintf(pFile, " : %E\n", m_Fwd[entryParameter]);
    }

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
SMP_Program()

Optimize the model using the sampling algorithm.
**************************************************************************************************************************************************************/
void SMP_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the sampling algorithm
       SamplingAlgorithm *Sampler = new SamplingAlgorithm();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       Sampler->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
} /* end SMP_Program() */
