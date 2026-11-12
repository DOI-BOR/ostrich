/**************************************************************************************************************************************************************

Powell's method is a zero-order (derivative-free) optimization algorithm, which utilizes the concept of conjugate directions to determine the optimal search
direction. Each iteration performs a line search along each of the current search directions, followed by a line search along the constructed conjugate
direction. Line searches are performed by the OptSearchClass2 class, with model evaluations solved using the shared, asynchronous MPI work-farming
functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "PowellAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
PowellAlgorithm::PowellAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read Powell information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Powell block within the file
        if (CheckToken(pFile, "BeginPowellAlg", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndPowellAlg", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginPowellAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndPowellAlg") == NULL) {

                if (strstr(line, "ConvergenceVal") != NULL) {
                    // Set the convergence value
                    sscanf(line, "%s %lf", tmp, &m_ConvVal);
                }
                else if (strstr(line, "MaxIterations") != NULL) {
                    // Set the maximum number of iterations
                    sscanf(line, "%s %d", tmp, &m_MaxIter);
                }

                // Get the next line
                line = GetNxtDataLine(pFile, inputFilename);
            }
        }
        else {
            // The Powell algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Initialize the search directions and the line search. The line search solves its parameter sets through this algorithm.
    m_NumDirs = m_pParamGroup->GetNumParams();
    m_SearchDirs = std::vector<std::vector<double>>(m_NumDirs, std::vector<double>(m_NumDirs, 0.00));
    m_pSearch = new OptSearchClass2(m_pParamGroup, [this](std::vector<std::vector<double>> samples) { return Evaluate(samples); });

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the algorithm and it's member variables
**************************************************************************************************************************************************************/
void PowellAlgorithm::Destroy(void) {

    m_SearchDirs.clear();
    delete m_pSearch;
    m_pSearch = NULL;

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void PowellAlgorithm::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
Evaluate()

Solves a set of parameter sets and returns their objectives. The best solution of the analysis is updated so that the base class can preserve the best model.
**************************************************************************************************************************************************************/
std::vector<double> PowellAlgorithm::Evaluate(std::vector<std::vector<double>> samples) {

    std::vector<double> objectives(samples.size(), INFINITY);
    ManageSingleObjectiveIterations(samples, 0, m_NumDirs, objectives);

    for (int entrySample = 0; entrySample < (int)samples.size(); entrySample++) {
        if (objectives[entrySample] < m_BestObjective) {
            m_BestObjective = objectives[entrySample];
            m_BestAlternative = samples[entrySample];
        }
    }

    return objectives;
}

/**************************************************************************************************************************************************************
ResetSearchDirections()

Resets the search directions to the unit vectors.
**************************************************************************************************************************************************************/
void PowellAlgorithm::ResetSearchDirections(void) {

    for (int i = 0; i < m_NumDirs; i++) {
        for (int j = 0; j < m_NumDirs; j++) {
            m_SearchDirs[i][j] = (i == j) ? 1.00 : 0.00;
        }
    }
}

/**************************************************************************************************************************************************************
LineMinimize()

Performs a line minimization along the direction <Scur>, starting from <X> (which has an objective of <Fcur>). If the line search itself finds a new minimum,
the solution is first moved to that location. The optimum step is then applied to the solution, moving half the distance to any violated bound, and the new
location is evaluated. The move is kept if it does not increase the objective, otherwise the solution is reverted.

On return, <Scur> holds the scaled step, and <X> and <Fcur> hold the current solution. Returns true if the step was kept.
**************************************************************************************************************************************************************/
bool PowellAlgorithm::LineMinimize(std::vector<double>& X, double& Fcur, std::vector<double>& Scur) {

    // Line search for the optimum step size
    double fmin = Fcur;
    std::vector<double> pmin = X;
    double xmin = m_pSearch->CalcStepSize(X, Fcur, Scur, fmin, pmin);

    if (fmin < Fcur) {
        // Found a new minimum when computing the step size
        X = pmin;
        Fcur = fmin;
    }

    // Apply the step
    std::vector<double> Xold = X;
    for (int j = 0; j < m_NumDirs; j++) {
        double upr = m_pParamGroup->GetParamPtr(j)->GetUpperBoundTransformed();
        double lwr = m_pParamGroup->GetParamPtr(j)->GetLowerBoundTransformed();

        Scur[j] *= xmin;
        X[j] += Scur[j];

        // If it's out of bounds, move half the distance to upr/lwr
        if (X[j] > upr) { X[j] = (upr + Xold[j]) / 2.00; m_NumUprViols++; }
        if (X[j] < lwr) { X[j] = (Xold[j] + lwr) / 2.00; m_NumLwrViols++; }
    }

    // Evaluate the new location
    double Ftmp = Evaluate(std::vector<std::vector<double>>{ X })[0];
    m_AlgCount++;

    if (Ftmp <= Fcur) {
        // Keep the move
        Fcur = Ftmp;
        return true;
    }

    // Revert if the move didn't improve the objective
    X = Xold;
    return false;
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using Powell's method.
**************************************************************************************************************************************************************/
void PowellAlgorithm::Optimize(void) {

    int n = m_NumDirs;
    std::vector<std::vector<double>>& S = m_SearchDirs;

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

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

    // Initialize the parameter vector
    std::vector<double> X;
    if (m_bWarmStart && (int)m_BestAlternative.size() >= n) {
        // Start from the best solution of the previous run
        X = std::vector<double>(m_BestAlternative.begin(), m_BestAlternative.begin() + n);
    }
    else {
        // Start from the initial parameter values
        for (int j = 0; j < n; j++) {
            X.push_back(m_pParamGroup->GetParamPtr(j)->GetInitialValueTransformed());
        }
    }

    if (!m_bWarmStart) {
        // Start a clean analysis
        // Write the information to the primary log file
        WriteSetup2("PowellAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), n, m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();
    }

    // Initialize the search directions
    ResetSearchDirections();

    // Evaluate the starting point
    double Fcur = Evaluate(std::vector<std::vector<double>>{ X })[0];
    double Fold = Fcur;
    m_AlgCount++;

    // Write the initial result to the log file
    SetParameterGroup(X);
    WriteRecord2(0, Fcur, Fcur, m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Main loop
    int count = 0;
    std::vector<double> Scur(n, 0.00);
    for (int iter = 0; iter < m_MaxIter; iter++) {

        m_CurIter = iter + 1;

        // Loop over all directions
        for (int i = 0; i < n; i++) {
            // Current direction
            Scur = S[i];

            // Normalize the search direction so that the largest value is 1.00, or reset the directions if the direction is empty
            double max = 0.00;
            for (int j = 0; j < n; j++) { max = MyMax(max, fabs(Scur[j])); }

            if (fabs(max) > NEARLY_ZERO) {
                for (int j = 0; j < n; j++) { Scur[j] /= max; }
            }
            else {
                ResetSearchDirections();
                Scur = S[i];
            }

            // Line minimization, saving the revised direction if the step was kept
            if (LineMinimize(X, Fcur, Scur)) {
                S[i] = Scur;
            }
        }

        // Construct the conjugate direction
        for (int j = 0; j < n; j++) {
            Scur[j] = 0.00;
            for (int i = 0; i < n; i++) { Scur[j] += S[i][j]; }
        }

        // Normalize the search direction so that the largest value is 1.00, if not zero
        double max = 0.00;
        for (int j = 0; j < n; j++) { max = MyMax(max, fabs(Scur[j])); }
        if (fabs(max) > NEARLY_ZERO) {
            for (int j = 0; j < n; j++) { Scur[j] /= max; }
        }

        // Line minimization along the conjugate direction
        if (LineMinimize(X, Fcur, Scur)) {
            // Shift the direction matrix and add the conjugate direction
            for (int i = 0; i < (n - 1); i++) { S[i] = S[i + 1]; }
            S[n - 1] = Scur;
        }
        else {
            // Restart the directions
            ResetSearchDirections();
            m_NumRestarts++;
        }

        // Write the iteration result to the log file
        SetParameterGroup(X);
        WriteRecord2(iter + 1, Fcur, fabs(Fold - Fcur), m_pObsGroup, m_pObjFunc, m_pParamGroup);

        // Convergence termination criteria
        if (fabs(Fold - Fcur) <= m_ConvVal) {
            count++;
            if (count >= m_MaxStallCount) { break; }
        }
        else {
            count = 0;
        }
        Fold = Fcur;
    }

    // Report the final solution as the optimal. Moves that tie the current objective are accepted by the algorithm, so the final solution can differ from
    // the first parameter set found with the same objective.
    if (Fcur <= m_BestObjective) {
        m_BestObjective = Fcur;
        m_BestAlternative = X;
    }

    // Terminate the secondary workers and let them clean up
    TerminateWorkers();

    // Write algorithm metrics
    WriteEndingMetrics();
}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void PowellAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Iterations          : %d\n", m_MaxIter);
    fprintf(pFile, "Convergence Val         : %E\n", m_ConvVal);

    // Write the parameter header to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "dObjFunc\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void PowellAlgorithm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    SetParameterGroup(m_BestAlternative);
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm        : Powell's Method (Conjugate Directions)\n");
    fprintf(pFile, "Max Iterations   : %d\n", m_MaxIter);
    fprintf(pFile, "Convergence Val  : %E\n", m_ConvVal);
    fprintf(pFile, "Iterations       : %d\n", m_CurIter);
    fprintf(pFile, "Algorithm Evals  : %d\n", m_AlgCount);
    fprintf(pFile, "Total Evals      : %d\n", m_NumSolves);
    fprintf(pFile, "Alg. Restarts    : %d\n", m_NumRestarts);
    fprintf(pFile, "Upper Violations : %d\n", m_NumUprViols);
    fprintf(pFile, "Lower Violations : %d\n", m_NumLwrViols);
    m_pSearch->WriteMetrics(pFile);

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
PWL_Program()

Optimize the model using Powell's method.
**************************************************************************************************************************************************************/
void PWL_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        // Create the Powell algorithm
        PowellAlgorithm* PowAlg = new PowellAlgorithm();

        // Begin the primary worker operations. This transfers information to the secondary workers
        // and manages the solution process
        PowAlg->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        // Create and setup the secondary worker, telling it to configure from MPI
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
    }
} /* end PWL_Program() */
