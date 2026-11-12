/**************************************************************************************************************************************************************

The Fletcher-Reeves algorithm is a first-order optimization algorithm, which utilizes the concept of conjugate directions in conjunction with the steepest-
descent information (negative of the gradient). Each iteration performs a line search along the current conjugate direction, computes the gradient at the
new location, and updates the conjugate direction using the Polak-Ribiere formula (reverting to steepest descent if the direction is not a descent
direction). The gradient is computed by the OptMathClass2 class and the line searches are performed by the OptSearchClass2 class, with model evaluations
solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "FletchReevesAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
FletchReevesAlgorithm::FletchReevesAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read Fletcher-Reeves information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Fletcher-Reeves block within the file
        if (CheckToken(pFile, "BeginFletchReevesAlg", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndFletchReevesAlg", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginFletchReevesAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndFletchReevesAlg") == NULL) {

                if (strstr(line, "ConvergenceVal") != NULL) {
                    // Set the convergence value
                    sscanf(line, "%s %lf", tmp, &m_ConvVal);
                }
                else if (strstr(line, "MaxStalls") != NULL) {
                    // Set the number of converged iterations required to stop
                    sscanf(line, "%s %d", tmp, &m_MaxCount);
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
            // The Fletcher-Reeves algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Initialize the gradient and line search classes. Both solve their parameter sets through this algorithm.
    m_NumParams = m_pParamGroup->GetNumParams();
    m_pMath = new OptMathClass2(m_pParamGroup, [this](std::vector<std::vector<double>> samples) { return Evaluate(samples); });
    m_pSearch = new OptSearchClass2(m_pParamGroup, [this](std::vector<std::vector<double>> samples) { return Evaluate(samples); });

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the algorithm and it's member variables
**************************************************************************************************************************************************************/
void FletchReevesAlgorithm::Destroy(void) {

    delete m_pMath;
    m_pMath = NULL;
    delete m_pSearch;
    m_pSearch = NULL;

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void FletchReevesAlgorithm::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
Evaluate()

Solves a set of parameter sets and returns their objectives. The best solution of the analysis is updated so that the base class can preserve the best model.
Empty parameter sets are not solved and are returned with an infinite objective.
**************************************************************************************************************************************************************/
std::vector<double> FletchReevesAlgorithm::Evaluate(std::vector<std::vector<double>> samples) {

    std::vector<double> objectives(samples.size(), INFINITY);
    ManageSingleObjectiveIterations(samples, 0, m_NumParams, objectives);

    for (int entrySample = 0; entrySample < (int)samples.size(); entrySample++) {
        if (!samples[entrySample].empty() && (objectives[entrySample] < m_BestObjective)) {
            m_BestObjective = objectives[entrySample];
            m_BestAlternative = samples[entrySample];
        }
    }

    return objectives;
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using the Fletcher-Reeves algorithm.
**************************************************************************************************************************************************************/
void FletchReevesAlgorithm::Optimize(void) {

    int n = m_NumParams;

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
    std::vector<double> p;
    if (m_bWarmStart && (int)m_BestAlternative.size() >= n) {
        // Start from the best solution of the previous run
        p = std::vector<double>(m_BestAlternative.begin(), m_BestAlternative.begin() + n);
    }
    else {
        // Start from the initial parameter values
        for (int j = 0; j < n; j++) {
            p.push_back(m_pParamGroup->GetParamPtr(j)->GetInitialValueTransformed());
        }
    }

    if (!m_bWarmStart) {
        // Start a clean analysis
        // Write the information to the primary log file
        WriteSetup2("FletchReevesAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), n, m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();
    }

    // Get the upper and lower bounds
    std::vector<double> upr, lwr;
    for (int j = 0; j < n; j++) {
        upr.push_back(m_pParamGroup->GetParamPtr(j)->GetUpperBoundTransformed());
        lwr.push_back(m_pParamGroup->GetParamPtr(j)->GetLowerBoundTransformed());
    }

    // Evaluate the starting point
    double fret = Evaluate(std::vector<std::vector<double>>{ p })[0];
    double old_fret = fret;
    double fmin = fret;
    std::vector<double> pmin = p;
    m_AlgCount++;

    // Write the initial result to the log file
    SetParameterGroup(p);
    WriteRecord2(0, fret, fret, m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Initialize the gradient
    std::vector<double> grad = m_pMath->CalcGradient(p, fret, fmin, pmin);
    if (fmin < fret) {
        // Found a better minimum during the gradient calculations
        p = pmin;
        fret = fmin;
        old_fret = fret;
    }

    // Initialize the search directions
    std::vector<double> g(n), h(n), xi(n);
    for (int j = 0; j < n; j++) {
        g[j] = -grad[j];
        xi[j] = h[j] = g[j];
    }

    int count = 0;
    for (int its = 0; its < m_MaxIter; its++) {

        // Normalize the search direction so that the largest value is 1.00. A zero direction means the starting gradient vanished, so the search is
        // complete before any iteration is performed.
        double max = 0.00;
        for (int j = 0; j < n; j++) { max = MyMax(max, fabs(xi[j])); }
        if (max <= 0.00) { break; }
        for (int j = 0; j < n; j++) { xi[j] /= max; }

        m_CurIter = its + 1;

        // Line minimization
        if (fret < fmin) { fmin = fret; pmin = p; }
        double xmin = m_pSearch->CalcStepSize(p, fret, xi, fmin, pmin);

        // Apply the step
        std::vector<double> old_p = p;
        for (int j = 0; j < n; j++) {
            xi[j] *= xmin;
            p[j] += xi[j];

            // If it's out of bounds, move half the distance to upr/lwr
            if (p[j] > upr[j]) { p[j] = (upr[j] + old_p[j]) / 2.00; m_NumUprViols++; }
            if (p[j] < lwr[j]) { p[j] = (old_p[j] + lwr[j]) / 2.00; m_NumLwrViols++; }
        }
        double ftmp = Evaluate(std::vector<std::vector<double>>{ p })[0];
        m_AlgCount++;

        if (ftmp <= fret) {
            fret = ftmp;
        }
        else {
            // Revert if the move didn't improve the objective
            p = old_p;
        }

        if (fmin < fret) {
            // The line search found a new optimal
            p = pmin;
            fret = fmin;
        }

        // Compute the gradient at the new location
        if (fret < fmin) { fmin = fret; pmin = p; }
        grad = m_pMath->CalcGradient(p, fret, fmin, pmin);
        if (fmin < fret) {
            // Found a better minimum during the gradient calculations
            p = pmin;
            fret = fmin;
        }

        // If the gradient vanished, the solution is at a stationary point and the search is complete. The direction is not updated, since a zero
        // gradient provides no search information (and would otherwise be counted as a steepest-descent restart).
        bool bZeroGradient = true;
        for (int j = 0; j < n; j++) {
            if (grad[j] != 0.00) { bZeroGradient = false; }
        }
        if (bZeroGradient) {
            SetParameterGroup(p);
            WriteRecord2(m_CurIter, fret, fabs(old_fret - fret), m_pObsGroup, m_pObjFunc, m_pParamGroup);
            break;
        }

        // Update the conjugate direction (Polak-Ribiere)
        double gg = 0.00;
        double dgg = 0.00;
        for (int j = 0; j < n; j++) {
            xi[j] = grad[j];
            gg += g[j] * g[j];
            dgg += (xi[j] + g[j]) * xi[j];
        }

        if (gg != 0.00) {
            double gam = dgg / gg;
            for (int j = 0; j < n; j++) {
                g[j] = -xi[j];
                xi[j] = h[j] = g[j] + gam * h[j];
            }

            // Check the slope, using steepest descent if the direction is not a descent direction
            double m = DotProduct(h.data(), grad.data(), n);
            if (m >= 0) {
                for (int j = 0; j < n; j++) { xi[j] = h[j] = -grad[j]; }
                m_NumRestarts++;
            }
        }

        // Write the iteration result to the log file
        SetParameterGroup(p);
        WriteRecord2(m_CurIter, fret, fabs(old_fret - fret), m_pObsGroup, m_pObjFunc, m_pParamGroup);

        // If the gradient (gg) is 0.00, we are done
        if (gg == 0.00) { break; }

        // Convergence termination criteria
        if (fabs(old_fret - fret) <= m_ConvVal) {
            count++;
            if (count >= m_MaxCount) { break; }
        }
        else {
            count = 0;
        }
        old_fret = fret;
    }

    // Report the final solution as the optimal. Moves that tie the current objective are accepted by the algorithm, so the final solution can differ from
    // the first parameter set found with the same objective.
    if (fret <= m_BestObjective) {
        m_BestObjective = fret;
        m_BestAlternative = p;
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
void FletchReevesAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Iterations          : %d\n", m_MaxIter);
    fprintf(pFile, "Max Stalls              : %d\n", m_MaxCount);
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
void FletchReevesAlgorithm::WriteEndingMetrics(void) {

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
    fprintf(pFile, "Algorithm        : Fletcher-Reeves (Conjugate Gradient)\n");
    fprintf(pFile, "Max Iterations   : %d\n", m_MaxIter);
    fprintf(pFile, "Convergence Val  : %E\n", m_ConvVal);
    fprintf(pFile, "Iterations       : %d\n", m_CurIter);
    fprintf(pFile, "Algorithm Evals  : %d\n", m_AlgCount);
    fprintf(pFile, "Total Evals      : %d\n", m_NumSolves);
    fprintf(pFile, "S-D Restarts     : %d\n", m_NumRestarts);
    fprintf(pFile, "Upper Violations : %d\n", m_NumUprViols);
    fprintf(pFile, "Lower Violations : %d\n", m_NumLwrViols);
    m_pMath->WriteMetrics(pFile);
    m_pSearch->WriteMetrics(pFile);

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
FLRV_Program()

Optimize the model using the Fletcher-Reeves algorithm.
**************************************************************************************************************************************************************/
void FLRV_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        // Create the Fletcher-Reeves algorithm
        FletchReevesAlgorithm* FRAlg = new FletchReevesAlgorithm();

        // Begin the primary worker operations. This transfers information to the secondary workers
        // and manages the solution process
        FRAlg->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        // Create and setup the secondary worker, telling it to configure from MPI
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
    }
} /* end FLRV_Program() */
