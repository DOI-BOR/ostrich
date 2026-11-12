/**************************************************************************************************************************************************************

Generalized Likelihood Uncertainty Estimation (GLUE), as described in Beven and Binley (1992). GLUE repeatedly draws batches of uniformly random parameter
sets and retains the best NumBehavioral samples. Samples with an objective below the behavioral threshold are considered behavioral. The search ends when
the desired number of behavioral samples has been found, or when the sample budget (MaxSamples) has been consumed. Each batch of samples is solved using
the MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "GLUE2.h"
#include "WriteUtility2.h"

/**************************************************************************************************************************************************************
CTOR

Reads the GLUE configuration from the BeginGLUE block of the OSTRICH input file. Default values are used for any settings that are not specified.
**************************************************************************************************************************************************************/
GLUE::GLUE() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read the GLUE configuration from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile == NULL) {
        // Couldn't open the file. Use the defaults and log the error.
        LogError(ERR_FILE_IO, "Couldn't open GLUE config. file. Using Defaults");
    }
    else {
        // Identify the GLUE block within the file
        if (CheckToken(pFile, "BeginGLUE", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndGLUE", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginGLUE", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndGLUE") == NULL) {
                if (strstr(line, "SamplesPerIter") != NULL) {
                    sscanf(line, "%s %d", tmp, &m_SamplesPerIter);
                    if (m_SamplesPerIter < 1) {
                        LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 10.");
                        m_SamplesPerIter = 10;
                    }
                }
                else if (strstr(line, "NumBehavioral") != NULL) {
                    sscanf(line, "%s %d", tmp, &m_NumDesired);
                    if (m_NumDesired < 1) {
                        LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 10.");
                        m_NumDesired = 10;
                    }
                }
                else if (strstr(line, "MaxSamples") != NULL) {
                    sscanf(line, "%s %lld", tmp, &m_MaxSamples);
                    if (m_MaxSamples < 1) {
                        LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 100.");
                        m_MaxSamples = 100;
                    }
                }
                else if (strstr(line, "Threshold") != NULL) {
                    sscanf(line, "%s %lf", tmp, &m_Threshold);
                }
                else {
                    sprintf(tmp, "Unknown token: %s", line);
                    LogError(ERR_FILE_IO, tmp);
                }

                // Get the next line
                line = GetNxtDataLine(pFile, inputFilename);
            }
        }
        else {
            // The GLUE algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default GLUE setup.");
        }

        fclose(pFile);
    }

    // Determine the maximum number of iterations. Whole iterations of SamplesPerIter samples are performed until the sample budget is consumed.
    m_MaxIters = 1 + (int)(m_MaxSamples / (long long)m_SamplesPerIter);

    // Initialize the list of retained samples
    m_Behavioral = std::vector<Sample>(m_NumDesired);

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the GLUE algorithm and it's member variables
**************************************************************************************************************************************************************/
void GLUE::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
Optimize()

Approximate the parameter uncertainty using GLUE. Batches of uniformly random samples are generated and solved, and the list of retained samples is
revised after each batch.
**************************************************************************************************************************************************************/
void GLUE::Optimize(void) {

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Warm start is not available for GLUE. Log the request and continue with a new analysis.
    if (m_bWarmStart) {
        LogError(ERR_FILE_IO, "Warm start is not available for the GLUE algorithm. A new analysis will be performed.");
        std::cout << "Warm start is not available for the GLUE algorithm. A new analysis will be performed." << std::endl;
    }

    // Fill the solve counter with zeros (no warm start)
    std::vector<int> solveCounts(numberOfMpiProcesses, 0);

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Write the information to the primary log file
    WriteSetup2("GLUE", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();

    // Main loop. Each iteration generates and solves a batch of SamplesPerIter samples.
    for (int entryIter = 0; entryIter < m_MaxIters; entryIter++) {
        // Check for termination
        if (IsQuit() == true) { break; }
        if (m_NumFound == m_NumDesired) { break; }
        long long numSamples = (long long)entryIter * (long long)m_SamplesPerIter;
        if (numSamples >= m_MaxSamples) { break; }

        m_CurIter = entryIter + 1;

        // Generate the batch of random samples
        std::vector<std::vector<double>> samples;
        for (int entrySample = 0; entrySample < m_SamplesPerIter; entrySample++) {
            samples.push_back(GenerateSample());
        }

        // Evaluate the batch. This is dispatched across all available MPI workers by the base class.
        std::vector<double> objectives(m_SamplesPerIter, INFINITY);
        ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);

        // Revise the list of retained samples and count the behavioral samples
        UpdateBehavioral(samples, objectives);

        // The first retained sample is always the best
        m_BestObjective = m_Behavioral[0].objective;
        m_BestAlternative = m_Behavioral[0].parameters;

        // Write the iteration result to the log file
        SetParameterGroup(m_BestAlternative);
        WriteRecord2(m_CurIter, m_BestObjective, (double)m_NumFound, m_pObsGroup, m_pObjFunc, m_pParamGroup);
    }

    // Terminate the secondary workers and let them clean up.
    TerminateWorkers();

    // Write algorithm metrics
    SetParameterGroup(m_BestAlternative);
    WriteEndingMetrics();

} /* end Optimize() */

/**************************************************************************************************************************************************************
GenerateSample()

Generate a parameter set by sampling each parameter uniformly within its (transformed) bounds. The sample is passed through the parameter group so that
integer parameters hold the integer values that will be evaluated.
**************************************************************************************************************************************************************/
std::vector<double> GLUE::GenerateSample(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    std::vector<double> sample;

    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
        double lwr = temp->GetLowerBoundTransformed();
        double upr = temp->GetUpperBoundTransformed();

        double r = (double)MyRand() / (double)MY_RAND_MAX;
        sample.push_back((r * (upr - lwr)) + lwr);
    }

    // Apply the sample to the parameter group and read back the values (integer parameters are rounded)
    SetParameterGroup(sample);
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        sample[entryParameter] = m_pParamGroup->GetParamPtr(entryParameter)->GetEstimatedValueTransformed();
    }

    return sample;
}

/**************************************************************************************************************************************************************
UpdateBehavioral()

Insert each solved sample into the sorted list of retained samples, if it improves upon a retained sample, and count the number of retained samples that
are behavioral. Samples are processed in the order they were generated; a sample must be strictly better than a retained sample to be inserted ahead of it.
**************************************************************************************************************************************************************/
void GLUE::UpdateBehavioral(std::vector<std::vector<double>>& samples, std::vector<double>& objectives) {

    for (size_t entrySample = 0; entrySample < samples.size(); entrySample++) {
        for (int entryRetained = 0; entryRetained < m_NumDesired; entryRetained++) {
            if (objectives[entrySample] < m_Behavioral[entryRetained].objective) {
                // Shift the remainder of the list to make room, discarding the last entry
                for (int entryShift = m_NumDesired - 1; entryShift > entryRetained; entryShift--) {
                    m_Behavioral[entryShift] = m_Behavioral[entryShift - 1];
                }

                // Insert the new entry
                m_Behavioral[entryRetained].objective = objectives[entrySample];
                m_Behavioral[entryRetained].parameters = samples[entrySample];
                break;
            }
        }
    }

    // Count the number of behavioral samples
    m_NumFound = 0;
    for (int entryRetained = 0; entryRetained < m_NumDesired; entryRetained++) {
        if (m_Behavioral[entryRetained].objective < m_Threshold) {
            m_NumFound++;
        }
    }
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Set a parameter set (transformed values) into the parameter group
**************************************************************************************************************************************************************/
void GLUE::SetParameterGroup(std::vector<double> parameters) {

    for (size_t entryParameter = 0; entryParameter < parameters.size(); entryParameter++) {
        m_pParamGroup->GetParamPtr((int)entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out the GLUE setup and the header of the run record
**************************************************************************************************************************************************************/
void GLUE::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Samples Per Iteration   : %d\n", m_SamplesPerIter);
    fprintf(pFile, "Desired Num. Behavioral : %d\n", m_NumDesired);
    fprintf(pFile, "Max Samples             : %lld\n", m_MaxSamples);
    fprintf(pFile, "Behavioral Threshold    : %E\n", m_Threshold);

    // Write the run record header to the file
    fprintf(pFile, "Gen   ");
    fprintf(pFile, "%-12s  ", "Best Value");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Num Found\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out the optimal solution, the algorithm metrics and the retained samples
**************************************************************************************************************************************************************/
void GLUE::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Generalized Likelihood Uncertainty Estimation\n");
    fprintf(pFile, "Behavioral Threshold    : %E\n", m_Threshold);
    fprintf(pFile, "Max Samples             : %lld\n", m_MaxSamples);
    fprintf(pFile, "Actual Num. Behavioral  : %d\n", m_NumFound);
    fprintf(pFile, "Desired Num. Behavioral : %d\n", m_NumDesired);
    fprintf(pFile, "Iterations              : %d\n", m_CurIter);
    fprintf(pFile, "Total Evals             : %d\n\n", m_NumSolves);

    // Write the retained samples, best first
    fprintf(pFile, "Sample  obj.function  ");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "\n");

    for (int entryRetained = 0; entryRetained < m_NumDesired; entryRetained++) {
        // Skip entries that have not been filled (fewer samples were solved than NumBehavioral)
        if (m_Behavioral[entryRetained].parameters.empty()) { continue; }

        fprintf(pFile, "%-4d  ", entryRetained);
        WritePreciseNumber2(pFile, m_Behavioral[entryRetained].objective);
        fprintf(pFile, "  ");
        SetParameterGroup(m_Behavioral[entryRetained].parameters);
        m_pParamGroup->Write(pFile, WRITE_SCI);
        fprintf(pFile, "\n");
    }

    // Restore the best solution into the parameter group
    SetParameterGroup(m_BestAlternative);

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
GLUE_Program()

Approximate the parameter uncertainty using GLUE.
**************************************************************************************************************************************************************/
void GLUE_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        GLUE* pGLUE = new GLUE();

        // Begin the primary worker operations. This transfers information to the secondary workers and manages the solution process
        pGLUE->Optimize();
    }
    else {
        // This is a secondary worker. Create the secondary worker, telling it to configure from MPI, and start the work
        ModelWorker worker = ModelWorker(true);
        worker.WorkMPI();
    }
} /* end GLUE_Program() */
