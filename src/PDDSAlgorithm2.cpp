/**************************************************************************************************************************************************************

Parallel Dynamically Dimensioned Search (PDDS) algorithm. PDDS is a parallel version of the DDS algorithm described in Tolson and Shoemaker (2007), Water
Resources Research.

The pre-refactor implementation farmed candidates out to the secondary workers asynchronously, generating each new candidate from the best solution known
when a worker became free. The Algorithm base class solves sets of parameter sets, so this version generates a batch of candidates from the current best
solution (one per available worker), solves the batch in parallel, and then processes the results in the order that the candidates were generated.

Original C++ port by L. Shawn Matott (2014) of Bryan Tolson's parallel Fortran DDS code.

**************************************************************************************************************************************************************/

// Include the class header
#include "PDDSAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values. The BeginParallelDDSAlg,
BeginParallelDDS, BeginParaDDSAlg, BeginParaDDS, BeginDDSAlg, and BeginDDS group headings are accepted, in that order of precedence.
**************************************************************************************************************************************************************/
PDDSAlgorithm::PDDSAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char beginToken[DEF_STR_SZ];
    char endToken[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read PDDS information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile == NULL) {
        FileOpenFailure("PDDSAlgorithm::CTOR", inputFilename);
    }

    // Accept multiple section headings
    const char* headings[] = { "ParallelDDSAlg", "ParallelDDS", "ParaDDSAlg", "ParaDDS", "DDSAlg", "DDS" };
    strcpy(beginToken, "BeginDDS");
    strcpy(endToken, "EndDDS");
    for (const char* heading : headings) {
        rewind(pFile);
        sprintf(tmp, "Begin%s", heading);
        if (CheckToken(pFile, tmp, inputFilename) == true) {
            sprintf(beginToken, "Begin%s", heading);
            sprintf(endToken, "End%s", heading);
            break;
        }
    }
    rewind(pFile);

    // Identify the PDDS block within the file
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
            else if (strstr(line, "UseOpt") != NULL) {
                // Set the special user option
                char option[DEF_STR_SZ];
                sscanf(line, "%s %s", tmp, option);
                m_UseOpt = option;
            }
            else if (strstr(line, "AlphaValue") != NULL) {
                // Set the perturbation probability used by the no-rand-num option
                sscanf(line, "%s %lf", tmp, &m_alpha);
            }
            else if (strstr(line, "BetaValue") != NULL) {
                // Set the perturbation multiplier used by the no-rand-num option
                sscanf(line, "%s %lf", tmp, &m_beta);
            }
            else if (strstr(line, "MaxIterations") != NULL) {
                // Set the maximum number of evaluations
                sscanf(line, "%s %d", tmp, &m_MaxIter);
            }
            else if (strstr(line, "UseInitialParamValues") != NULL) {
                // Start from the user-supplied parameter values
                m_UserSuppliedInit = true;
            }
            else if (strstr(line, "UseRandomParamValues") != NULL) {
                // Start from random parameter values
                m_UserSuppliedInit = false;
            }
            else if (strstr(line, "EnableDebugging") != NULL) {
                // Write debugging information to the console
                m_bDebug = true;
            }

            // Get the next line
            line = GetNxtDataLine(pFile, inputFilename);
        }
    }
    else {
        // The PDDS algorithm has been specified but there is no input block. Log an error.
        LogError(ERR_FILE_IO, "Using default PDDS algorithm setup.");
    }

    // Read any user-supplied initial parameter sets
    rewind(pFile);
    ReadInitParams(pFile, inputFilename);

    fclose(pFile);

    // Validate the configuration
    if ((m_r_val < 0) || (m_r_val > 1)) {
        LogError(ERR_FILE_IO, "Bad Perturbation value specified for PDDS Algorithm");
        ExitProgram(1);
    }
    if (m_MaxIter < 1) {
        LogError(ERR_FILE_IO, "Maximum PDDS Algorithm iterations must be >0");
        ExitProgram(1);
    }
    if ((m_UseOpt != "standard") && (m_UseOpt != "no-rand-num") && (m_UseOpt != "try-int-solution")) {
        LogError(ERR_FILE_IO, "Unknown UseOpt value specified for PDDS Algorithm. Using the standard option.");
        m_UseOpt = "standard";
    }

    // TODO: port support for the TelescopingStrategy option. Telescoping corrections are applied to each candidate by the pre-refactor algorithm,
    // but telescoping has not yet been re-enabled within the Algorithm base class.

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the algorithm and it's member variables
**************************************************************************************************************************************************************/
void PDDSAlgorithm::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
ReadInitParams()

Reads one or more initial parameter sets from the optional BeginInitParams/EndInitParams block. Each line contains a value for every parameter, in the order
the parameters are listed in the parameter group. Values are converted into the transformed space used by the search.
**************************************************************************************************************************************************************/
void PDDSAlgorithm::ReadInitParams(FILE* pFile, IroncladString inputFilename) {

    char* line;
    char tmp[DEF_STR_SZ];

    if (CheckToken(pFile, "BeginInitParams", inputFilename) == false) {
        return;
    }

    // Identify the end of the block
    FindToken(pFile, "EndInitParams", inputFilename);
    rewind(pFile);

    // Read each of the entries
    int numberOfParameters = m_pParamGroup->GetNumParams();
    FindToken(pFile, "BeginInitParams", inputFilename);
    line = GetNxtDataLine(pFile, inputFilename);

    while (strstr(line, "EndInitParams") == NULL) {
        char* pTok = line;
        std::vector<double> parameters;

        // Extract values, one-by-one, making any necessary conversions
        for (int k = 0; k < numberOfParameters; k++) {
            int j = ExtractString(pTok, tmp);
            j = ValidateExtraction(j, k, numberOfParameters, "PDDSAlgorithm::ReadInitParams()");
            pTok += j;
            parameters.push_back(m_pParamGroup->GetParamPtr(k)->ConvertInVal(atof(tmp)));
        }
        m_InitParams.push_back(parameters);

        line = GetNxtDataLine(pFile, inputFilename);
    }
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void PDDSAlgorithm::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
CorrectParameters()

Applies any user-specified parameter corrections (expert judgement) to a parameter set and returns the corrected values.
**************************************************************************************************************************************************************/
std::vector<double> PDDSAlgorithm::CorrectParameters(std::vector<double> parameters) {

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
std::vector<double> PDDSAlgorithm::Evaluate(std::vector<std::vector<double>> samples) {

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

Generates a neighboring decision variable value for a single decision variable being perturbed by the PDDS algorithm. With the standard option the value is
perturbed using a normal distribution with a standard deviation of r times the parameter range. The no-rand-num option scales the current value by the beta
value, and the try-int-solution option truncates the current value to an integer. The new value respects the parameter bounds, using reflection at the bounds.

Ported to C from neigh_value.f90 (Bryan Tolson, Nov 2005).
**************************************************************************************************************************************************************/
double PDDSAlgorithm::PerturbParam(double x_cur, ParameterABC* pParam) {

    double x_new;
    double x_max = pParam->GetUpperBoundTransformed();
    double x_min = pParam->GetLowerBoundTransformed();

    if (m_UseOpt == "no-rand-num") {
        x_new = x_cur * m_beta;
    }
    else if (m_UseOpt == "try-int-solution") {
        x_new = (double)((int)(x_cur));
    }
    else {
        x_new = x_cur + GaussRandom() * m_r_val * (x_max - x_min);
    }

    // If the new value is outside of the bounds, reflect it back inside
    if (x_new < x_min) {
        x_new = x_min + (x_min - x_new);

        // If the reflection goes past x_max then the value should be x_min, since without reflection the approach goes well past the lower bound. This
        // keeps x close to the lower bound when x_cur is close to the lower bound.
        if (x_new > x_max) { x_new = x_min; }
    }
    else if (x_new > x_max) {
        x_new = x_max - (x_new - x_max);

        if (x_new < x_min) { x_new = x_max; }
    }

    return x_new;
}

/**************************************************************************************************************************************************************
WriteStatus()

Appends the current best solution to the dds_status.out file so that the last best result is available during long runs.
**************************************************************************************************************************************************************/
void PDDSAlgorithm::WriteStatus(int evaluation, double objective, std::vector<double> parameters) {

    FILE* pOut = fopen("dds_status.out", "a");
    if (pOut == NULL) { return; }

    fprintf(pOut, "%4d\t%E", evaluation, objective);
    for (int k = 0; k < (int)parameters.size(); k++) {
        fprintf(pOut, "\t%E", parameters[k]);
    }
    fprintf(pOut, "\n");
    fclose(pOut);
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using PDDS.
**************************************************************************************************************************************************************/
void PDDSAlgorithm::Optimize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Determine the number of candidates that can be solved concurrently
    m_BatchSize = (numberOfMpiProcesses - 1) + (m_bSolveOnPrimary ? 1 : 0);
    if (m_BatchSize < 1) {
        m_BatchSize = 1;
    }

    // Perform initial warm start recovery, if needed
    std::vector<int> solveCounts;

    if (m_bWarmStart) {
        // Read worker files for the previous solves
        std::vector<std::vector<std::vector<double>>> workerValues = ReadWorkerFiles(numberOfMpiProcesses, m_bSolveOnPrimary);

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

    if (!m_bWarmStart) {
        // Start a clean analysis
        // Write the information to the primary log file
        WriteSetup2("PDDSAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();
    }

    // Initialize the file that records the probability of perturbing each decision variable
    FILE* pPnFile = fopen("OstDDSPn.txt", "w");
    fprintf(pPnFile, "EVAL  Pn\n");
    fclose(pPnFile);

    // Initialize the status file that records each improvement to the best solution
    FILE* pStatusFile = fopen("dds_status.out", "w");
    fprintf(pStatusFile, "STEP\tOBJ._FUNCTION");
    for (int k = 0; k < numberOfParameters; k++) {
        fprintf(pStatusFile, "\t%-12s", m_pParamGroup->GetParamPtr(k)->GetName());
    }
    fprintf(pStatusFile, "\n");
    fclose(pStatusFile);

    //------------------------------------------------------------
    //				INITIALIZATION
    //------------------------------------------------------------

    // Determine the starting point(s) and the number of evaluations used for initialization. At least one evaluation per worker is used, and any
    // evaluations that are not covered by a starting point are filled with uniform random samples.
    bool bWarmStart = (m_bWarmStart && ((int)m_BestAlternative.size() >= numberOfParameters));
    std::vector<std::vector<double>> startingPoints;
    int initFunctEvals;

    if (bWarmStart) {
        // Start from the best solution of the previous run
        startingPoints.push_back(std::vector<double>(m_BestAlternative.begin(), m_BestAlternative.begin() + numberOfParameters));
        initFunctEvals = iMax(m_BatchSize, 1);
    }
    else if (m_UserSuppliedInit) {
        // Start from the user-supplied initial solutions. Use the initial parameter values if no BeginInitParams block was provided.
        startingPoints = m_InitParams;
        if (startingPoints.empty()) {
            std::vector<double> initialValues;
            for (int k = 0; k < numberOfParameters; k++) {
                initialValues.push_back(m_pParamGroup->GetParamPtr(k)->GetInitialValueTransformed());
            }
            startingPoints.push_back(initialValues);
        }
        initFunctEvals = iMax(m_BatchSize, (int)startingPoints.size());
    }
    else {
        // Standard random initialization as described in Tolson and Shoemaker (2007). Use the best of at least 5 (or one per worker) uniform random samples.
        initFunctEvals = iMax(iMax(5, (int)(0.005 * double(m_MaxIter))), m_BatchSize);
    }

    // The initialization cannot exceed the budget
    if (initFunctEvals > m_MaxIter) {
        initFunctEvals = m_MaxIter;
    }

    if (m_bDebug) {
        printf("ini_fevals, ileft = %d , %d\n", initFunctEvals, m_MaxIter);
    }

    // Create the initialization samples and solve them together
    std::vector<std::vector<double>> samples;
    for (int i = 0; i < initFunctEvals; i++) {
        std::vector<double> testParams;

        if (i < (int)startingPoints.size()) {
            testParams = startingPoints[i];
        }
        else {
            for (int k = 0; k < numberOfParameters; k++) {
                ParameterABC* pParam = m_pParamGroup->GetParamPtr(k);
                testParams.push_back(pParam->GetLowerBoundTransformed() +
                                     UniformRandom() * (pParam->GetUpperBoundTransformed() - pParam->GetLowerBoundTransformed()));
            }
        }

        samples.push_back(CorrectParameters(testParams));
    }

    std::vector<double> objectives = Evaluate(samples);

    // Select the best initialization sample, processing the samples in order
    double Fbest = INFINITY;
    std::vector<double> bestParams;
    bool bLastRecorded = false;
    for (int i = 1; i <= initFunctEvals; i++) {
        double Ftest = objectives[i - 1];

        if ((i == 1) || (Ftest <= Fbest)) {
            // Initializing, or the solution is the best. Update the solution.
            Fbest = Ftest;
            bestParams = samples[i - 1];

            if (m_bDebug) { printf("%4d\t%E\n", i, Fbest); }
            WriteStatus(i, Fbest, bestParams);
            SetParameterGroup(bestParams);
            WriteRecord2(i, Fbest, (double)(m_MaxIter - i), m_pObsGroup, m_pObjFunc, m_pParamGroup);
            bLastRecorded = true;
        }
        else {
            bLastRecorded = false;
        }
    }

    //------------------------------------------------------------
    //				MAIN PDDS LOOP
    //------------------------------------------------------------

    // The first evaluation of each worker after the initialization perturbs every decision variable. The perturbation probability then decreases
    // logarithmically over the remaining budget.
    int numberOfFullPerturbations = 2 * m_BatchSize;
    double logRemaining = log(double(m_MaxIter - numberOfFullPerturbations));

    int evaluation = initFunctEvals + 1;
    while (evaluation <= m_MaxIter) {

        // Generate a batch of candidates from the current best solution
        int batchSize = std::min(m_BatchSize, m_MaxIter - evaluation + 1);
        std::vector<std::vector<double>> candidates;

        for (int b = 0; b < batchSize; b++) {
            int eval = evaluation + b;

            // Determine the probability that a decision variable is selected for perturbation
            double Pn;
            if (m_UseOpt == "no-rand-num") {
                Pn = m_alpha;
            }
            else if (eval <= numberOfFullPerturbations) {
                Pn = 1.00;
            }
            else if (logRemaining > 0.00) {
                Pn = 1.0 - log(double(eval - numberOfFullPerturbations)) / logRemaining;
            }
            else {
                Pn = 0.00;
            }

            if (m_bDebug) { printf("Pn = %e\n", Pn); }
            pPnFile = fopen("OstDDSPn.txt", "a");
            fprintf(pPnFile, "%04d  %E\n", eval, Pn);
            fclose(pPnFile);

            // Define the test parameters initially as the current best solution
            std::vector<double> testParams = bestParams;

            // Randomly select decision variables to perturb
            int dvnCount = 0;
            for (int k = 0; k < numberOfParameters; k++) {
                if (UniformRandom() < Pn) {
                    dvnCount++;
                    testParams[k] = PerturbParam(bestParams[k], m_pParamGroup->GetParamPtr(k));
                }
            }

            // If no decision variables were selected at random, select one
            if (dvnCount == 0) {
                int dv = (int)(ceil((double)(numberOfParameters) * UniformRandom())) - 1;
                dv = std::max(0, std::min(dv, numberOfParameters - 1));
                testParams[dv] = PerturbParam(bestParams[dv], m_pParamGroup->GetParamPtr(dv));
            }

            candidates.push_back(CorrectParameters(testParams));
        }

        // Solve the batch of candidates
        objectives = Evaluate(candidates);

        // Process the results in the order that the candidates were generated
        for (int b = 0; b < batchSize; b++) {
            int eval = evaluation + b;
            double Ftest = objectives[b];

            if (Ftest <= Fbest) {
                // Update the current (best) solution
                Fbest = Ftest;
                bestParams = candidates[b];

                if (m_bDebug) { printf("%4d\t%E\n", eval, Fbest); }
                WriteStatus(eval, Fbest, bestParams);
                SetParameterGroup(bestParams);
                WriteRecord2(eval, Fbest, (double)(m_MaxIter - eval), m_pObsGroup, m_pObjFunc, m_pParamGroup);
                bLastRecorded = true;
            }
            else {
                bLastRecorded = false;
            }
        }

        evaluation += batchSize;
    }

    // Record the final state of the search, if it was not already recorded by the last evaluation
    if (!bLastRecorded) {
        SetParameterGroup(bestParams);
        WriteRecord2(m_MaxIter, Fbest, 0.00, m_pObsGroup, m_pObjFunc, m_pParamGroup);
    }

    // Report the final solution as the optimal. Candidates that tie the current best are accepted by PDDS, so the final solution can differ from the
    // first parameter set found with the same objective.
    if (Fbest <= m_BestObjective) {
        m_BestObjective = Fbest;
        m_BestAlternative = bestParams;
    }

    // Terminate the secondary workers and let them clean up
    TerminateWorkers();

    // Write algorithm metrics
    WriteEndingMetrics();

    if (m_bDebug) { printf("DDS optimization is complete!\n"); }
}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void PDDSAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Iterations          : %d\n", m_MaxIter);
    fprintf(pFile, "Perturbation Value      : %lf\n", m_r_val);
    fprintf(pFile, "Batch Size              : %d\n", m_BatchSize);
    fprintf(pFile, "Initialization          : ");
    if (m_bWarmStart) { fprintf(pFile, "Warm Start\n"); }
    else if (m_UserSuppliedInit) { fprintf(pFile, "User Supplied (%d parameter sets)\n", iMax(1, (int)m_InitParams.size())); }
    else { fprintf(pFile, "Random Parameter Values\n"); }

    // Write the parameter header to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Trials Remaining\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void PDDSAlgorithm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Write the result information
    SetParameterGroup(m_BestAlternative);
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Parallel Dynamically-Dimensioned Search Algorithm (PDDS)\n");
    fprintf(pFile, "Desired Convergence Val : N/A\n");
    fprintf(pFile, "Actual Convergence Val  : N/A\n");
    fprintf(pFile, "Max Generations         : %d\n", m_MaxIter);
    fprintf(pFile, "Actual Generations      : %d\n", m_CurIter);
    fprintf(pFile, "Peterbation Value       : %lf\n", m_r_val);
    fprintf(pFile, "Debug Statements        : %s\n", m_bDebug ? "enabled" : "disabled");
    fprintf(pFile, "Initial Solution        : ");
    if (m_bWarmStart) { fprintf(pFile, "Warm Start\n"); }
    else if (m_UserSuppliedInit) { fprintf(pFile, "User Supplied\n"); }
    else { fprintf(pFile, "Randomly Generated\n"); }
    fprintf(pFile, "Special User Option     : %s\n", m_UseOpt.c_str());
    if (m_UseOpt == "no-rand-num") {
        fprintf(pFile, "Alpha Value             : %lf\n", m_alpha);
        fprintf(pFile, "Beta Value              : %lf\n", m_beta);
    }
    else {
        fprintf(pFile, "Alpha Value             : not used\n");
        fprintf(pFile, "Beta Value              : not used\n");
    }
    fprintf(pFile, "Number of Processors    : %d\n", numberOfMpiProcesses);
    fprintf(pFile, "Batch Size              : %d\n", m_BatchSize);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);
    fprintf(pFile, "Algorithm successfully converged on a solution, however more runs may be needed\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
PDDS_Program()

Optimize the model using PDDS.
**************************************************************************************************************************************************************/
void PDDS_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        // Create the PDDS algorithm
        PDDSAlgorithm* PDDS = new PDDSAlgorithm();

        // Begin the primary worker operations. This transfers information to the secondary workers
        // and manages the solution process
        PDDS->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        // Create and setup the secondary worker, telling it to configure from MPI
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
    }
} /* end PDDS_Program() */
