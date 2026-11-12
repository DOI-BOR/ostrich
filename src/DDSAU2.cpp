/**************************************************************************************************************************************************************

DDS for Approximation of Uncertainty (DDS-AU), as described in Tolson and Shoemaker (2008), Water Resources Research. DDS-AU performs a number of independent,
randomly initialized, short-budget DDS searches and selects a single parameter set from each search to approximate the parameter uncertainty.

The pre-refactor implementation ran each search with a separate DDS (or parallel DDS) algorithm object. The Algorithm base class configures and terminates the
MPI workers once per analysis, so this version performs the DDS searches within a single algorithm. When ParallelSearches is enabled, the searches are
advanced together, with one candidate from every active search solved concurrently at each step. Otherwise, the searches are performed one after another.

Original implementation by L. Shawn Matott, based on the DDS-AU method of Bryan Tolson.

**************************************************************************************************************************************************************/

// Include the class header
#include "DDSAU2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values within the
Begin_DDSAU_Alg group.
**************************************************************************************************************************************************************/
DDSAU::DDSAU() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char yesno[DEF_STR_SZ];
    const char* beginToken = "Begin_DDSAU_Alg";
    const char* endToken = "End_DDSAU_Alg";
    IroncladString inputFilename = GetInFileName();

    // Read DDS-AU information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile == NULL) {
        FileOpenFailure("DDSAU::CTOR", inputFilename);
    }

    // Identify the DDS-AU block within the file
    if (CheckToken(pFile, beginToken, inputFilename) == true) {
        // Identify the end of the algorithm parameter block
        FindToken(pFile, endToken, inputFilename);
        rewind(pFile);

        // Import data starting from the beginning of the block
        FindToken(pFile, beginToken, inputFilename);
        line = GetNxtDataLine(pFile, inputFilename);

        while (strstr(line, endToken) == NULL) {

            if (strstr(line, "PerturbationValue") != NULL) {
                // Set the perturbation value
                sscanf(line, "%s %lf", tmp, &m_r_val);
            }
            else if (strstr(line, "NumSearches") != NULL) {
                // Set the number of independent searches
                sscanf(line, "%s %d", tmp, &m_nsols);
            }
            else if (strstr(line, "Threshold") != NULL) {
                // Set the behavioral threshold
                sscanf(line, "%s %lf", tmp, &m_fmax);
            }
            else if (strstr(line, "MinItersPerSearch") != NULL) {
                // Set the minimum budget of each search
                sscanf(line, "%s %d", tmp, &m_MinIter);
            }
            else if (strstr(line, "MaxItersPerSearch") != NULL) {
                // Set the maximum budget of each search
                sscanf(line, "%s %d", tmp, &m_MaxIter);
            }
            else if (strstr(line, "ParallelSearches") != NULL) {
                // Advance the searches together
                sscanf(line, "%s %s", tmp, yesno);
                MyStrLwr(yesno);
                m_bParallel = (strcmp(yesno, "yes") == 0);
            }
            else if (strstr(line, "Randomize") != NULL) {
                // Randomly select a behavioral solution from each search
                sscanf(line, "%s %s", tmp, yesno);
                MyStrLwr(yesno);
                m_bRandomize = (strcmp(yesno, "yes") == 0);
            }
            else if (strstr(line, "ReviseAU") != NULL) {
                // Re-use the results of previously completed searches
                sscanf(line, "%s %s", tmp, yesno);
                MyStrLwr(yesno);
                m_bReviseAU = (strcmp(yesno, "yes") == 0);
            }

            // Get the next line
            line = GetNxtDataLine(pFile, inputFilename);
        }
    }
    else {
        // The DDS-AU algorithm has been specified but there is no input block. Log an error.
        LogError(ERR_FILE_IO, "Using default DDSAU algorithm setup.");
    }

    fclose(pFile);

    // Validate the configuration
    if ((m_r_val < 0) || (m_r_val > 1)) {
        LogError(ERR_FILE_IO, "Bad Perturbation value specified for DDSAU Algorithm");
        ExitProgram(1);
    }
    if (m_nsols < 1) {
        LogError(ERR_FILE_IO, "Number of DDSAU searches must be >0");
        ExitProgram(1);
    }
    if (m_MaxIter < 1) {
        LogError(ERR_FILE_IO, "Maximum DDSAU Algorithm iterations must be >0");
        ExitProgram(1);
    }
    if (m_MinIter < 1) {
        LogError(ERR_FILE_IO, "Minimum DDSAU Algorithm iterations must be >0");
        ExitProgram(1);
    }
    if (m_MaxIter < m_MinIter) {
        int itmp = m_MaxIter;
        m_MaxIter = m_MinIter;
        m_MinIter = itmp;
    }

    // TODO: port support for the TelescopingStrategy option. Telescoping corrections are applied to each DDS candidate by the pre-refactor algorithm,
    // but telescoping has not yet been re-enabled within the Algorithm base class.

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the algorithm and it's member variables
**************************************************************************************************************************************************************/
void DDSAU::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void DDSAU::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
CorrectParameters()

Applies any user-specified parameter corrections (expert judgement) to a parameter set and returns the corrected values.
**************************************************************************************************************************************************************/
std::vector<double> DDSAU::CorrectParameters(std::vector<double> parameters) {

    SetParameterGroup(parameters);
    PerformParameterCorrections();

    std::vector<double> corrected;
    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        corrected.push_back(m_pParamGroup->GetParamPtr(entryParameter)->GetEstimatedValueTransformed());
    }

    return corrected;
}

/**************************************************************************************************************************************************************
Evaluate()

Solves a set of parameter sets and returns their objectives. The best solution of the analysis is updated so that the base class can preserve the best model.
**************************************************************************************************************************************************************/
std::vector<double> DDSAU::Evaluate(std::vector<std::vector<double>> samples) {

    std::vector<double> objectives(samples.size(), INFINITY);
    ManageSingleObjectiveIterations(samples, 0, m_pParamGroup->GetNumParams(), objectives);
    m_CurIter += (int)samples.size();

    for (int entrySample = 0; entrySample < (int)samples.size(); entrySample++) {
        if (objectives[entrySample] < m_BestObjective) {
            m_BestObjective = objectives[entrySample];
            m_BestAlternative = samples[entrySample];
        }
    }

    return objectives;
}

/**************************************************************************************************************************************************************
PerturbParam()

Generates a neighboring decision variable value for a single decision variable, following the DDS algorithm. Real-valued parameters are perturbed using a
normal distribution with a standard deviation of r times the parameter range. Integer parameters are perturbed using a uniform distribution and are
guaranteed to change by at least one unit. The new value respects the parameter bounds, using reflection at the bounds.
**************************************************************************************************************************************************************/
double DDSAU::PerturbParam(double x_best, ParameterABC* pParam) {

    double x_new, r;
    double x_max = pParam->GetUpperBoundTransformed();
    double x_min = pParam->GetLowerBoundTransformed();

    if (strcmp(pParam->GetType(), "integer") == 0) {
        r = 2.00 * UniformRandom() - 1.00; // -1 < r < +1
        x_new = x_best + r * m_r_val * (x_max - x_min);

        // Check that a discrete change has been made
        if ((int)x_new == (int)x_best) {
            if (x_new < x_best) x_new = x_new - 1.00; // decrement
            else x_new = x_new + 1.00;                // increment
        }

        // Bounds check
        if ((int)x_new > (int)x_max) x_new = x_best - 1.00;
        if ((int)x_new < (int)x_min) x_new = x_best + 1.00;
    }
    else {
        // Real-valued (and all other) parameters
        x_new = x_best + GaussRandom() * m_r_val * (x_max - x_min);
    }

    // If the new value is outside of the bounds, reflect it back inside
    if (x_new < x_min) {
        x_new = x_min + (x_min - x_new);
        if (x_new > x_max) { x_new = x_min; }
    }
    else if (x_new > x_max) {
        x_new = x_max - (x_new - x_max);
        if (x_new < x_min) { x_new = x_max; }
    }

    return x_new;
}

/**************************************************************************************************************************************************************
RunSearches()

Performs a set of independent DDS searches together. Each search is randomly initialized and then follows the DDS algorithm for its own budget. At each step,
one candidate from every active search is generated and the candidates are solved concurrently. The results are then processed search by search, so each
search is unaffected by the others.
**************************************************************************************************************************************************************/
void DDSAU::RunSearches(std::vector<int> searchIndices) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Processes the result of an evaluation for a search, accepting the candidate if it is no worse than the current best
    auto processResult = [&](Search& search, std::vector<double> parameters, double objective) {
        search.objectives.push_back(objective);
        search.parameters.push_back(parameters);
        int evaluation = (int)search.objectives.size();

        if ((evaluation == 1) || (objective <= search.Fbest)) {
            search.Fbest = objective;
            search.bestParams = parameters;
            search.recordEvals.push_back(evaluation);
        }
    };

    //------------------------------------------------------------
    //				INITIALIZATION
    //------------------------------------------------------------

    // Create the random initialization samples of every search and solve them together
    std::vector<std::vector<double>> samples;
    std::vector<int> owners;
    for (int s : searchIndices) {
        for (int i = 0; i < m_Searches[s].initFunctEvals; i++) {
            std::vector<double> testParams;
            for (int k = 0; k < numberOfParameters; k++) {
                ParameterABC* pParam = m_pParamGroup->GetParamPtr(k);
                testParams.push_back(UniformRandom() * (pParam->GetUpperBoundTransformed() - pParam->GetLowerBoundTransformed()) +
                                     pParam->GetLowerBoundTransformed());
            }
            samples.push_back(CorrectParameters(testParams));
            owners.push_back(s);
        }
    }

    std::vector<double> objectives = Evaluate(samples);
    for (int entry = 0; entry < (int)samples.size(); entry++) {
        processResult(m_Searches[owners[entry]], samples[entry], objectives[entry]);
    }

    //------------------------------------------------------------
    //				MAIN DDS LOOP
    //------------------------------------------------------------

    // Determine the longest search
    int maxRemaining = 0;
    for (int s : searchIndices) {
        maxRemaining = std::max(maxRemaining, m_Searches[s].budget - m_Searches[s].initFunctEvals);
    }

    for (int t = 1; t <= maxRemaining; t++) {

        // Generate a candidate for every search that has budget remaining
        samples.clear();
        owners.clear();

        for (int s : searchIndices) {
            Search& search = m_Searches[s];
            int itersRemaining = search.budget - search.initFunctEvals;
            if (t > itersRemaining) {
                continue;
            }

            // Determine the probability that a decision variable is selected as a neighbor
            double Pn = (itersRemaining > 1) ? (1.0 - log(double(t)) / log(double(itersRemaining))) : 0.00;

            // Define the test parameters initially as the current best solution
            std::vector<double> testParams = search.bestParams;

            // Randomly select decision variables in the neighborhood to perturb
            int dvnCount = 0;
            for (int k = 0; k < numberOfParameters; k++) {
                if (UniformRandom() < Pn) {
                    dvnCount++;
                    testParams[k] = PerturbParam(search.bestParams[k], m_pParamGroup->GetParamPtr(k));
                }
            }

            // If no decision variables were selected at random, select one
            if (dvnCount == 0) {
                int dv = (int)(ceil((double)(numberOfParameters) * UniformRandom())) - 1;
                dv = std::max(0, std::min(dv, numberOfParameters - 1));
                testParams[dv] = PerturbParam(search.bestParams[dv], m_pParamGroup->GetParamPtr(dv));
            }

            samples.push_back(CorrectParameters(testParams));
            owners.push_back(s);
        }

        // Solve the candidates and update each search
        objectives = Evaluate(samples);
        for (int entry = 0; entry < (int)samples.size(); entry++) {
            processResult(m_Searches[owners[entry]], samples[entry], objectives[entry]);
        }
    }
}

/**************************************************************************************************************************************************************
SelectSolution()

Selects the evaluation that represents a search. The best solution of the search is used unless Randomize is enabled and the search found at least one
behavioral solution, in which case one of the behavioral solutions is selected at random.
**************************************************************************************************************************************************************/
void DDSAU::SelectSolution(Search& search) {

    // Count the behavioral solutions and identify the best solution
    int nBehavioral = 0;
    int bestIndex = 0;
    for (int entry = 0; entry < (int)search.objectives.size(); entry++) {
        if (search.objectives[entry] <= m_fmax) { nBehavioral++; }
        if (search.objectives[entry] < search.objectives[bestIndex]) { bestIndex = entry; }
    }

    if ((m_bRandomize == false) || (nBehavioral == 0)) {
        // Use the best solution
        search.selected = bestIndex;
    }
    else {
        // Randomly select a behavioral solution
        int whichBehavioral = 1 + (MyRand() % nBehavioral);
        int count = 0;
        for (int entry = 0; entry < (int)search.objectives.size(); entry++) {
            if (search.objectives[entry] <= m_fmax) { count++; }
            if (count == whichBehavioral) {
                search.selected = entry;
                break;
            }
        }
    }
}

/**************************************************************************************************************************************************************
ReadSearch()

Reads the evaluations of a previously completed search from the OstModel0_DDS[M].txt file. Returns false if the file is missing or contains no evaluations,
in which case the search is performed again.
**************************************************************************************************************************************************************/
bool DDSAU::ReadSearch(int searchIndex) {

    char fileName[DEF_STR_SZ];
    sprintf(fileName, "OstModel0_DDS%d.txt", searchIndex);
    FILE* pFile = fopen(fileName, "r");
    if (pFile == NULL) {
        return false;
    }

    int numberOfParameters = m_pParamGroup->GetNumParams();
    Search search;
    char line[DEF_STR_SZ * 4];

    // Skip the header and read each evaluation
    bool bHeader = true;
    while (fgets(line, sizeof(line), pFile) != NULL) {
        if (bHeader) { bHeader = false; continue; }

        std::vector<double> values;
        char* pTok = strtok(line, " \t\r\n");
        while (pTok != NULL) {
            values.push_back(atof(pTok));
            pTok = strtok(NULL, " \t\r\n");
        }
        if ((int)values.size() < numberOfParameters + 2) {
            continue;
        }

        // Convert the parameters into the transformed space used by the search
        std::vector<double> parameters;
        for (int k = 0; k < numberOfParameters; k++) {
            parameters.push_back(m_pParamGroup->GetParamPtr(k)->ConvertInVal(values[k + 2]));
        }

        search.objectives.push_back(values[1]);
        search.parameters.push_back(parameters);
        int evaluation = (int)search.objectives.size();
        if ((evaluation == 1) || (values[1] <= search.Fbest)) {
            search.Fbest = values[1];
            search.bestParams = parameters;
            search.recordEvals.push_back(evaluation);
        }
    }
    fclose(pFile);

    if (search.objectives.empty()) {
        return false;
    }

    printf("Using previous results located in %s\n", fileName);
    search.budget = (int)search.objectives.size();
    search.bRevised = true;
    m_Searches[searchIndex] = search;

    // Track the best solution of the analysis
    for (int entry = 0; entry < (int)search.objectives.size(); entry++) {
        if (search.objectives[entry] < m_BestObjective) {
            m_BestObjective = search.objectives[entry];
            m_BestAlternative = search.parameters[entry];
        }
    }

    return true;
}

/**************************************************************************************************************************************************************
Optimize()

Approximate the parameter uncertainty using DDS-AU.
**************************************************************************************************************************************************************/
void DDSAU::Optimize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Warm start is not used by DDS-AU. Previously completed searches can be re-used with the ReviseAU option.
    if (m_bWarmStart) {
        LogError(ERR_FILE_IO, "WarmStart is not supported by DDSAU and will be ignored. Use the ReviseAU option to re-use previous searches.");
        m_bWarmStart = false;
    }

    // Initialize the workers
    ConfigureWorkers(std::vector<int>(numberOfMpiProcesses, 0));

    // Write the information to the primary log file
    WriteSetup2("DDSAU", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();

    // Delete the results of previous searches, unless they are being revised
    if (!m_bReviseAU) {
        char fileName[DEF_STR_SZ];
        for (int proc = 0; proc < numberOfMpiProcesses; proc++) {
            for (int i = 0; i < m_nsols; i++) {
                sprintf(fileName, "OstModel%d_DDS%d.txt", proc, i);
                remove(fileName);
                sprintf(fileName, "OstOutput%d_DDS%d.txt", proc, i);
                remove(fileName);
            }
        }
    }

    // Configure each search
    m_Searches = std::vector<Search>(m_nsols);
    std::vector<int> searchesToRun;
    for (int i = 0; i < m_nsols; i++) {
        // Determine the budget of the search
        Search& search = m_Searches[i];
        if (m_MinIter == m_MaxIter) {
            search.budget = m_MaxIter;
        }
        else {
            search.budget = m_MinIter + (MyRand() % ((m_MaxIter - m_MinIter) + 1));
        }

        // Determine the number of random initialization samples, as in the DDS algorithm
        search.initFunctEvals = iMax(5, (int)(0.005 * double(search.budget)));
        if (search.initFunctEvals >= search.budget) {
            search.initFunctEvals = 1;
        }

        // Re-use the results of a previous search, if requested and available
        if (m_bReviseAU && ReadSearch(i)) {
            continue;
        }
        searchesToRun.push_back(i);
    }

    // Perform the searches
    if (m_bParallel) {
        RunSearches(searchesToRun);
    }
    else {
        for (int i : searchesToRun) {
            RunSearches(std::vector<int>{ i });
        }
    }

    // Terminate the secondary workers and let them clean up
    TerminateWorkers();

    // Post-process the searches
    m_nbhvr = 0;
    for (int i = 0; i < m_nsols; i++) {
        Search& search = m_Searches[i];
        SelectSolution(search);
        if (search.objectives[search.selected] <= m_fmax) {
            m_nbhvr++;
        }

        if (!search.bRevised) {
            WriteSearch(i);
        }

        // Record the solution selected from the search
        SetParameterGroup(search.parameters[search.selected]);
        WriteRecord2(i, search.objectives[search.selected], (double)search.budget, m_pObsGroup, m_pObjFunc, m_pParamGroup);
    }

    // Write the summary and algorithm metrics
    WriteSummary();
    WriteEndingMetrics();
}

/**************************************************************************************************************************************************************
WriteSolutionLine()

Writes an evaluation number, objective, and parameter set (in the untransformed space) to a file.
**************************************************************************************************************************************************************/
void DDSAU::WriteSolutionLine(FILE* pFile, int run, double objective, std::vector<double> parameters) {

    fprintf(pFile, "%-4d  %E  ", run, objective);
    SetParameterGroup(parameters);
    m_pParamGroup->Write(pFile, WRITE_SCI);
    fprintf(pFile, "\n");
}

/**************************************************************************************************************************************************************
WriteSearch()

Writes the output files of a search. OstModel0_DDS[M].txt lists every evaluation of the search and is read by the ReviseAU option. OstOutput0_DDS[M].txt lists
each improvement to the best solution of the search, followed by a final record.
**************************************************************************************************************************************************************/
void DDSAU::WriteSearch(int searchIndex) {

    Search& search = m_Searches[searchIndex];
    char fileName[DEF_STR_SZ];

    // Write every evaluation of the search
    sprintf(fileName, "OstModel0_DDS%d.txt", searchIndex);
    FILE* pFile = fopen(fileName, "w");
    fprintf(pFile, "Run   obj.function  ");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "\n");
    for (int entry = 0; entry < (int)search.objectives.size(); entry++) {
        WriteSolutionLine(pFile, entry + 1, search.objectives[entry], search.parameters[entry]);
    }
    fclose(pFile);

    // Write the improvements to the best solution of the search
    sprintf(fileName, "OstOutput0_DDS%d.txt", searchIndex);
    pFile = fopen(fileName, "w");
    fprintf(pFile, "DDS Search %d of the DDS-AU analysis\n", searchIndex);
    fprintf(pFile, "Budget                  : %d\n", search.budget);
    fprintf(pFile, "Initialization Samples  : %d\n", search.initFunctEvals);
    fprintf(pFile, "Run   %-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Trials Remaining\n");
    for (int evaluation : search.recordEvals) {
        fprintf(pFile, "%-4d  %E  ", evaluation, search.objectives[evaluation - 1]);
        SetParameterGroup(search.parameters[evaluation - 1]);
        m_pParamGroup->Write(pFile, WRITE_SCI);
        fprintf(pFile, "%E\n", (double)(search.budget - evaluation));
    }
    if (search.recordEvals.empty() || (search.recordEvals.back() != search.budget)) {
        fprintf(pFile, "%-4d  %E  ", search.budget, search.Fbest);
        SetParameterGroup(search.bestParams);
        m_pParamGroup->Write(pFile, WRITE_SCI);
        fprintf(pFile, "%E\n", 0.00);
    }
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteSummary()

Writes the OstOutputDDSAU.txt summary file, which lists the solution selected from each search, the algorithm metrics, and the behavioral solutions.
**************************************************************************************************************************************************************/
void DDSAU::WriteSummary(void) {

    FILE* pFile = fopen("OstOutputDDSAU.txt", "w");
    if (m_bParallel) { fprintf(pFile, "Parallel DDS for Approximation of Uncertainty (PDDSAU)\n"); }
    else { fprintf(pFile, "DDS for Approximation of Uncertainty (DDSAU)\n"); }

    // Write the solution selected from each search
    fprintf(pFile, "Iter  Run   obj.function  ");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "\n");
    for (int i = 0; i < m_nsols; i++) {
        Search& search = m_Searches[i];
        fprintf(pFile, "%-4d  ", i);
        WriteSolutionLine(pFile, search.selected + 1, search.objectives[search.selected], search.parameters[search.selected]);
    }

    // Write the algorithm metrics
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm                : DDS for Approximating Uncertainty (DDSAU)\n");
    fprintf(pFile, "Perturbation Value       : %0.2lf\n", m_r_val);
    fprintf(pFile, "Desired # of Samples     : %d\n", m_nsols);
    fprintf(pFile, "Actual # of Samples      : %d\n", m_nbhvr);
    fprintf(pFile, "Min DDS Evals per Sample : %d\n", m_MinIter);
    fprintf(pFile, "Max DDS Evals per Sample : %d\n", m_MaxIter);
    fprintf(pFile, "Behavioral Threshold     : %E\n", m_fmax);
    fprintf(pFile, "Parallel searches?       : %s\n", m_bParallel ? "yes" : "no");
    fprintf(pFile, "Randomize samples?       : %s\n", m_bRandomize ? "yes" : "no");
    fprintf(pFile, "Revise Previous DDS AU?  : %s\n", m_bReviseAU ? "yes" : "no");
    fprintf(pFile, "Total Evals              : %d\n", m_NumSolves);

    // Write the behavioral solutions
    fprintf(pFile, "\nList of Behavioral Solutions\n");
    fprintf(pFile, "Iter  Run   obj.function  ");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "\n");
    for (int i = 0; i < m_nsols; i++) {
        Search& search = m_Searches[i];
        if (search.objectives[search.selected] <= m_fmax) {
            fprintf(pFile, "%-4d  ", i);
            WriteSolutionLine(pFile, search.selected + 1, search.objectives[search.selected], search.parameters[search.selected]);
        }
    }

    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void DDSAU::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Number of Searches      : %d\n", m_nsols);
    fprintf(pFile, "Evals per Search        : %d to %d\n", m_MinIter, m_MaxIter);
    fprintf(pFile, "Perturbation Value      : %lf\n", m_r_val);
    fprintf(pFile, "Behavioral Threshold    : %E\n", m_fmax);
    fprintf(pFile, "Parallel Searches       : %s\n", m_bParallel ? "yes" : "no");
    fprintf(pFile, "Randomize Samples       : %s\n", m_bRandomize ? "yes" : "no");
    fprintf(pFile, "Revise Previous Results : %s\n", m_bReviseAU ? "yes" : "no");

    // Write the parameter header to the file. One record is written for the solution selected from each search.
    fprintf(pFile, "Search  ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Search Evals\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void DDSAU::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the best solution found by any search
    if ((int)m_BestAlternative.size() >= m_pParamGroup->GetNumParams()) {
        SetParameterGroup(m_BestAlternative);
    }
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : DDS for Approximating Uncertainty (DDSAU)\n");
    fprintf(pFile, "Desired # of Samples    : %d\n", m_nsols);
    fprintf(pFile, "Actual # of Samples     : %d\n", m_nbhvr);
    fprintf(pFile, "Behavioral Threshold    : %E\n", m_fmax);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);
    fprintf(pFile, "The selected and behavioral solutions are listed in OstOutputDDSAU.txt\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
DDSAU_Program()

Approximate the parameter uncertainty using DDS-AU.
**************************************************************************************************************************************************************/
void DDSAU_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        // Create the DDS-AU algorithm
        DDSAU* pDDSAU = new DDSAU();

        // Begin the primary worker operations. This transfers information to the secondary workers
        // and manages the solution process
        pDDSAU->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        // Create and setup the secondary worker, telling it to configure from MPI
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
    }
} /* end DDSAU_Program() */
