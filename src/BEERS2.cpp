/**************************************************************************************************************************************************************

Balanced Exploration-Exploitation Random Search (BEERS) applies a balanced approach to randomly search a parameter space. The search initially favors
exploration of the search space and will try to maximize the distance between points that are evaluated. As the search progresses the algorithm favors
exploitation and will favor points that are close to the current optimal. To facilitate exploration vs. exploitation the algorithm maintains an archive of
every point that is ever evaluated during the search.

Each candidate depends on the objective of every previously evaluated candidate, so the samples are evaluated one at a time using the shared MPI
work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "BEERS2.h"

/**************************************************************************************************************************************************************
CTOR

Reads the BEERS configuration from the input file. If the block is not present, the default number of samples is used.
**************************************************************************************************************************************************************/
BEERS::BEERS() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read the BEERS information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile == NULL) {
        LogError(ERR_FILE_IO, "Couldn't open BEERS config. file. Using Defaults");
    }
    else {
        // Identify the BEERS block within the file
        if (CheckToken(pFile, "BeginBEERS", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndBEERS", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginBEERS", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndBEERS") == NULL) {
                if (strstr(line, "NumSamples") != NULL) {
                    sscanf(line, "%s %d", tmp, &m_NumSamples);
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
            // The algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Validate the number of samples
    if (m_NumSamples < 1) {
        LogError(ERR_FILE_IO, "Invalid number of BEERS samples, defaulting to 1");
        m_NumSamples = 1;
    }

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by BEERS and it's member variables
**************************************************************************************************************************************************************/
void BEERS::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using BEERS.
**************************************************************************************************************************************************************/
void BEERS::Calibrate(void) {
    Optimize();
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using BEERS. A random first sample is evaluated, then each subsequent candidate is generated randomly and accepted with a
probability that blends exploration (distance to the nearest archived point) and exploitation (quality of the nearest archived point). The weight given
to exploitation increases linearly with the number of samples.
**************************************************************************************************************************************************************/
void BEERS::Optimize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Warm start is not available for BEERS. Log the request and continue with a new analysis.
    if (m_bWarmStart) {
        const char* msg = "Warm start is not available for the BEERS algorithm. A new analysis will be performed.";
        LogError(ERR_FILE_IO, msg);
        std::cout << msg << std::endl;
    }

    // Extracting the starting values from the model output is not available for BEERS
    if (m_pParamGroup->CheckExtraction()) {
        const char* msg = "Parameter extraction is not available for the BEERS algorithm. A random first sample will be used.";
        LogError(ERR_FILE_IO, msg);
        std::cout << msg << std::endl;
    }

    // Fill the solve counter with zeros
    std::vector<int> solveCounts = std::vector<int>(numberOfMpiProcesses, 0);

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    // Write the information to the primary log file
    WriteSetup2("BEERS", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();

    // Generate the random first sample and store the search space information
    m_Min = std::vector<double>(numberOfParameters, 0.00);
    m_Range = std::vector<double>(numberOfParameters, 0.00);

    ArchiveMember first;
    first.X = std::vector<double>(numberOfParameters, 0.00);
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
        double lwr = temp->GetLowerBoundTransformed();
        double upr = temp->GetUpperBoundTransformed();
        double range = upr - lwr;
        double r = (double)MyRand() / (double)MY_RAND_MAX;

        first.X[entryParameter] = (r * range) + lwr;
        m_Min[entryParameter] = lwr;
        m_Range[entryParameter] = range;
    }

    // Evaluate the first sample
    first.X = ApplyParameterGroup(first.X);
    first.F = Evaluate(first.X);
    m_Archive.push_back(first);
    m_BestIndex = 0;
    m_BestObjective = first.F;
    m_BestAlternative = first.X;
    m_CurSample = 1;

    // Write the initial configuration
    SetParameterGroup(m_BestAlternative);
    WriteRecord2(0, m_BestObjective, (double)(m_NumSamples - 1), m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Main optimization loop
    m_MaxDist = sqrt((double)numberOfParameters);
    for (int g = 1; g < m_NumSamples; g++) {

        // Assign model probabilities
        AssignModelProbs(m_Archive[m_BestIndex].F);

        WriteArchive();

        // Compute exploration and exploitation weights
        double Wexploit = ((double)(g + 1) / (double)m_NumSamples);
        double Wexplore = 1.00 - Wexploit;

        // Generate random candidates until one is accepted
        ArchiveMember candidate;
        candidate.X = std::vector<double>(numberOfParameters, 0.00);

        bool bCandidateAccepted = false;
        m_MinProbAccept = 0.00;
        while (bCandidateAccepted == false) {
            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                candidate.X[entryParameter] = (((double)MyRand() / (double)MY_RAND_MAX) * m_Range[entryParameter]) + m_Min[entryParameter];
            }

            // Compute exploitation and exploration probabilities
            double Pexploit, Pexplore;
            CalcProbabilities(candidate.X, Pexploit, Pexplore);

            double Paccept = (Wexploit * Pexploit) + (Wexplore * Pexplore);
            if (Paccept < m_MinProbAccept) {
                Paccept = m_MinProbAccept;
            }

            double Raccept = ((double)MyRand() / (double)MY_RAND_MAX);
            if (Raccept < Paccept) {
                bCandidateAccepted = true;
            }

            // Ensures that a candidate is eventually selected
            m_MinProbAccept += 1E-6;
        }

        // Evaluate the accepted candidate and add it to the archive
        candidate.X = ApplyParameterGroup(candidate.X);
        candidate.F = Evaluate(candidate.X);
        m_Archive.push_back(candidate);
        m_CurSample = g + 1;
        int currentIndex = (int)m_Archive.size() - 1;

        // Update the best solution
        if (candidate.F < m_Archive[m_BestIndex].F) {
            m_BestIndex = currentIndex;
            m_BestObjective = candidate.F;
            m_BestAlternative = candidate.X;

            SetParameterGroup(m_BestAlternative);
            WriteRecord2(g + 1, m_BestObjective, (double)(m_NumSamples - g - 1), m_pObsGroup, m_pObjFunc, m_pParamGroup);
        }

        AdjustRanks(currentIndex);
    }

    // Write the final archive, including the last evaluated sample
    AssignModelProbs(m_Archive[m_BestIndex].F);
    WriteArchive();

    // Terminate the secondary workers and let them clean up.
    TerminateWorkers();

    // Set the best solution into the parameter group and write the algorithm metrics
    SetParameterGroup(m_BestAlternative);
    WriteEndingMetrics();

} /* end Optimize() */

/**************************************************************************************************************************************************************
Evaluate()

Evaluates a single sample using the shared MPI work-farming functionality.
**************************************************************************************************************************************************************/
double BEERS::Evaluate(std::vector<double> parameters) {

    std::vector<std::vector<double>> samples = { parameters };
    std::vector<double> objectives(1, INFINITY);
    ManageSingleObjectiveIterations(samples, 0, (int)parameters.size(), objectives);

    return objectives[0];
}

/**************************************************************************************************************************************************************
ApplyParameterGroup()

Set a parameter set (transformed values) into the parameter group and read back the values. Integer parameters are rounded, so that the archived values
are the values that will be evaluated.
**************************************************************************************************************************************************************/
std::vector<double> BEERS::ApplyParameterGroup(std::vector<double> parameters) {

    SetParameterGroup(parameters);
    for (size_t entryParameter = 0; entryParameter < parameters.size(); entryParameter++) {
        parameters[entryParameter] = m_pParamGroup->GetParamPtr((int)entryParameter)->GetEstimatedValueTransformed();
    }
    return parameters;
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Set a parameter set (transformed values) into the parameter group
**************************************************************************************************************************************************************/
void BEERS::SetParameterGroup(std::vector<double> parameters) {

    for (size_t entryParameter = 0; entryParameter < parameters.size(); entryParameter++) {
        m_pParamGroup->GetParamPtr((int)entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
AdjustRanks()

Iterate through the archive and adjust ranks based on newest insertion.
**************************************************************************************************************************************************************/
void BEERS::AdjustRanks(int currentIndex) {

    ArchiveMember& current = m_Archive[currentIndex];
    for (int entryMember = 0; entryMember < (int)m_Archive.size(); entryMember++) {
        if (entryMember == currentIndex) continue;

        if (m_Archive[entryMember].F <= current.F) {
            current.Z++;
        }
        else {
            m_Archive[entryMember].Z++;
        }
    }
}

/**************************************************************************************************************************************************************
AssignModelProbs()

Compute the likelihood of each model in the current archive. If the best objective is not positive (or not finite, which happens when every evaluation
has failed), the probability is based on the rank instead of the objective ratio.
**************************************************************************************************************************************************************/
void BEERS::AssignModelProbs(double Fbest) {

    for (ArchiveMember& member : m_Archive) {
        if ((Fbest <= 0.00) || !std::isfinite(Fbest)) {
            member.P = 1.00 / (member.Z + 1);
        }
        else {
            member.P = Fbest / member.F;
        }
    }
}

/**************************************************************************************************************************************************************
CalcProbabilities()

Calculate the probability of accepting a candidate solution based on explorative and exploitative behavior.
**************************************************************************************************************************************************************/
void BEERS::CalcProbabilities(const std::vector<double>& candidate, double& Pexploit, double& Pexplore) {

    Pexplore = 1.00;
    Pexploit = 1.00;
    if (m_CurSample < 2) return;

    // Compute distance to nearest point in archive
    double dmin = NEARLY_HUGE;
    int numberOfParameters = (int)candidate.size();

    for (const ArchiveMember& member : m_Archive) {
        // Compute normalized distance between points
        double dtst = 0.00;
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            double x1i = (member.X[entryParameter] - m_Min[entryParameter]) / m_Range[entryParameter];
            double x2i = (candidate[entryParameter] - m_Min[entryParameter]) / m_Range[entryParameter];

            dtst += ((x2i - x1i) * (x2i - x1i));
        }
        dtst = sqrt(dtst);

        // Re-assign min distance, if needed
        if (dtst < dmin) {
            dmin = dtst;
            Pexploit = member.P;
        }
    }

    // Exploration acceptance probability is ratio of min and max distance
    Pexplore = (dmin / m_MaxDist);
}

/**************************************************************************************************************************************************************
WriteArchive()

Write out archive.
**************************************************************************************************************************************************************/
void BEERS::WriteArchive(void) {

    FILE* pOut = fopen("OstArchive.txt", "w");
    if (pOut == NULL) return;

    for (const ArchiveMember& member : m_Archive) {
        fprintf(pOut, "%E  ", member.F);
        fprintf(pOut, "%04d  ", (int)(member.Z));
        fprintf(pOut, "%E  ", member.P);
        for (size_t entryParameter = 0; entryParameter < member.X.size(); entryParameter++) {
            fprintf(pOut, "%E  ", member.X[entryParameter]);
        }
        fprintf(pOut, "\n");
    }
    fclose(pOut);
}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for BEERS.
**************************************************************************************************************************************************************/
void BEERS::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Desired Samples         : %d\n", m_NumSamples);

    // Write the parameter header to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Samples_Remaining\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void BEERS::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Balanced Exploration-Exploitation Random Search\n");
    fprintf(pFile, "Desired Samples         : %d\n", m_NumSamples);
    fprintf(pFile, "Actual Samples          : %d\n", m_CurSample);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
BEERS_Program()

Calibrate or optimize the model using BEERS.
**************************************************************************************************************************************************************/
void BEERS_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        BEERS* TheAlg = new BEERS();

        if (TheAlg->GetObjectiveFunctionType() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
        else { TheAlg->Optimize(); }

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor
        worker.WorkMPI();
    }
} /* end BEERS_Program() */
