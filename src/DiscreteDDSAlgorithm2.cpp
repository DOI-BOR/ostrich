/**************************************************************************************************************************************************************

Discrete Dynamically Dimensioned Search (DDDS) algorithm. A variant of the DDS algorithm of Tolson and Shoemaker (2007), Water Resources Research, intended for
discrete-valued decision variables. Each perturbation is forced to change the value of a decision variable by at least one unit.

The random initialization samples are independent and are solved together using the shared, asynchronous MPI work-farming functionality provided by the
Algorithm base class. The main search loop is inherently sequential, as each candidate is generated from the current best solution.

Original implementation by L. Shawn Matott (2011).

**************************************************************************************************************************************************************/

// Include the class header
#include "DiscreteDDSAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values. Values are read from
the BeginDiscreteDDSAlg group.
**************************************************************************************************************************************************************/
DiscreteDDSAlgorithm::DiscreteDDSAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char beginToken[DEF_STR_SZ];
    char endToken[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read DDS information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile == NULL) {
        FileOpenFailure("DiscreteDDSAlgorithm::CTOR", inputFilename);
    }

    // Set the section headings
    strcpy(beginToken, "BeginDiscreteDDSAlg");
    strcpy(endToken, "EndDiscreteDDSAlg");

    // Identify the Discrete DDS block within the file
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
            else if (strstr(line, "MaxIterations") != NULL) {
                // Set the maximum number of evaluations
                sscanf(line, "%s %d", tmp, &m_MaxIter);
            }
            else if (strstr(line, "UseInitialParamValues") != NULL) {
                // Start from the initial parameter values
                m_UserSuppliedInit = true;
            }
            else if (strstr(line, "UseRandomParamValues") != NULL) {
                // Start from random parameter values
                m_UserSuppliedInit = false;
            }

            // Get the next line
            line = GetNxtDataLine(pFile, inputFilename);
        }
    }
    else {
        // The Discrete DDS algorithm has been specified but there is no input block. Log an error.
        LogError(ERR_FILE_IO, "Using default Discrete DDS algorithm setup.");
    }

    fclose(pFile);

    // Validate the configuration
    if ((m_r_val < 0) || (m_r_val > 1)) {
        LogError(ERR_FILE_IO, "Bad Perturbation value specified for Discrete DDS Algorithm");
        ExitProgram(1);
    }
    if (m_MaxIter < 1) {
        LogError(ERR_FILE_IO, "Maximum Discrete DDS Algorithm iterations must be >0");
        ExitProgram(1);
    }

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the algorithm and it's member variables
**************************************************************************************************************************************************************/
void DiscreteDDSAlgorithm::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void DiscreteDDSAlgorithm::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
ApplyParameterGroup()

Sets a parameter set into the parameter group and returns the values as stored by the group. Integer parameters are rounded to the nearest integer, so the
values that are solved, recorded, and perturbed are the same. The pre-refactor algorithm achieved this by writing each candidate into the parameter group
before the model was executed.
**************************************************************************************************************************************************************/
std::vector<double> DiscreteDDSAlgorithm::ApplyParameterGroup(std::vector<double> parameters) {

    SetParameterGroup(parameters);

    std::vector<double> applied;
    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        applied.push_back(m_pParamGroup->GetParamPtr(entryParameter)->GetEstimatedValueTransformed());
    }

    return applied;
}

/**************************************************************************************************************************************************************
Evaluate()

Solves a set of parameter sets and returns their objectives. The best solution of the analysis is updated so that the base class can preserve the best model.
**************************************************************************************************************************************************************/
std::vector<double> DiscreteDDSAlgorithm::Evaluate(std::vector<std::vector<double>> samples) {

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

Generates a neighboring decision variable value for a single decision variable being perturbed by the Discrete DDS algorithm. The value is perturbed using a
normal distribution with a standard deviation of r times the parameter range. Since the decision variables are assumed to be discrete-valued, the new value
must differ from the current best by at least one unit. The new value respects the parameter bounds, using reflection at the bounds.
**************************************************************************************************************************************************************/
double DiscreteDDSAlgorithm::PerturbParam(double x_best, ParameterABC* pParam) {

    double x_max = pParam->GetUpperBoundTransformed();
    double x_min = pParam->GetLowerBoundTransformed();

    double x_new = x_best + GaussRandom() * m_r_val * (x_max - x_min);

    // Check that a discrete change has been made. Values are compared after rounding to the nearest integer, consistent with the rounding applied by
    // the parameter group, so that the shift below always changes the value by one unit.
    if (floor(x_new + 0.5) == floor(x_best + 0.5)) {
        if (x_new < x_best) x_new = x_new - 1.00; // decrement
        else x_new = x_new + 1.00;                // increment
        m_nCorr++;
    }

    // If the new value is outside of the bounds, reflect it back inside
    if (x_new < x_min) {
        x_new = x_min + (x_min - x_new);

        // If the reflection goes past x_max then the value should be x_min, since without reflection the approach goes well past the lower bound. This
        // keeps x close to the lower bound when x_best is close to the lower bound.
        if (x_new > x_max) { x_new = x_min; }
    }
    else if (x_new > x_max) {
        x_new = x_max - (x_new - x_max);

        if (x_new < x_min) { x_new = x_max; }
    }

    return x_new;
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using Discrete DDS.
**************************************************************************************************************************************************************/
void DiscreteDDSAlgorithm::Optimize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

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

    if (!m_bWarmStart) {
        // Start a clean analysis
        // Write the information to the primary log file
        WriteSetup2("DiscreteDDSAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();
    }

    //------------------------------------------------------------
    //				INITIALIZATION
    //------------------------------------------------------------

    // Determine the starting point(s) and the number of evaluations used for initialization
    bool bWarmStart = (m_bWarmStart && ((int)m_BestAlternative.size() >= numberOfParameters));
    std::vector<double> startingPoint;
    int initFunctEvals;

    if (bWarmStart) {
        // Start from the best solution of the previous run
        initFunctEvals = 1;
        startingPoint = std::vector<double>(m_BestAlternative.begin(), m_BestAlternative.begin() + numberOfParameters);
    }
    else if (m_UserSuppliedInit) {
        // Start from the user-supplied initial solution
        initFunctEvals = 1;
        for (int k = 0; k < numberOfParameters; k++) {
            startingPoint.push_back(m_pParamGroup->GetParamPtr(k)->GetInitialValueTransformed());
        }
    }
    else {
        // Standard random initialization. Calculate the number of evaluations to use for random sampling.
        initFunctEvals = iMax(5, (int)(0.005 * double(m_MaxIter)));
        if (initFunctEvals >= m_MaxIter) {
            initFunctEvals = 1;
        }
    }

    // Reduce the number of evaluations available to the main DDS loop
    int itersRemaining = m_MaxIter - initFunctEvals;
    if (itersRemaining <= 0) {
        LogError(ERR_FILE_IO, "DiscreteDDSAlgorithm: # of Initialization samples >= Max # func evaluations");
        ExitProgram(1);
    }

    // Create the initialization samples. Random samples are drawn in the same order as the pre-refactor algorithm, and are then solved together.
    std::vector<std::vector<double>> samples;
    for (int i = 1; i <= initFunctEvals; i++) {
        std::vector<double> testParams = startingPoint;

        if (!m_UserSuppliedInit && !bWarmStart) {
            testParams.clear();
            for (int k = 0; k < numberOfParameters; k++) {
                ParameterABC* pParam = m_pParamGroup->GetParamPtr(k);
                testParams.push_back(UniformRandom() * (pParam->GetUpperBoundTransformed() - pParam->GetLowerBoundTransformed()) +
                                     pParam->GetLowerBoundTransformed());
            }
        }

        samples.push_back(ApplyParameterGroup(testParams));
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

            SetParameterGroup(bestParams);
            WriteRecord2(i, Fbest, (double)(m_MaxIter - i), m_pObsGroup, m_pObjFunc, m_pParamGroup);
            bLastRecorded = (i == initFunctEvals);
        }
        else {
            bLastRecorded = false;
        }
    }

    //------------------------------------------------------------
    //				MAIN DISCRETE DDS LOOP
    //------------------------------------------------------------

    // Initialize the file that records the probability of perturbing each decision variable
    FILE* pPnFile = fopen("OstDDSPn.txt", "w");
    fprintf(pPnFile, "EVAL  Pn\n");
    fclose(pPnFile);

    for (int i = 1; i <= itersRemaining; i++) {

        // Determine the probability that a decision variable is selected as a neighbor
        double Pn = 1.0 - log(double(i)) / log(double(itersRemaining));
        pPnFile = fopen("OstDDSPn.txt", "a");
        fprintf(pPnFile, "%04d  %E\n", i, Pn);
        fclose(pPnFile);

        // Define the test parameters initially as the current best solution
        std::vector<double> testParams = bestParams;

        // Randomly select decision variables in the neighborhood to perturb
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

        // Evaluate the objective function for the candidate
        testParams = ApplyParameterGroup(testParams);
        double Ftest = Evaluate(std::vector<std::vector<double>>{ testParams })[0];

        if (Ftest <= Fbest) {
            // Update the current (best) solution
            Fbest = Ftest;
            bestParams = testParams;

            SetParameterGroup(bestParams);
            WriteRecord2(i + initFunctEvals, Fbest, (double)(m_MaxIter - (i + initFunctEvals)), m_pObsGroup, m_pObjFunc, m_pParamGroup);
            bLastRecorded = true;
        }
        else {
            bLastRecorded = false;
        }
    }

    // Record the final state of the search, if it was not already recorded by the last evaluation
    if (!bLastRecorded) {
        SetParameterGroup(bestParams);
        WriteRecord2(m_MaxIter, Fbest, 0.00, m_pObsGroup, m_pObjFunc, m_pParamGroup);
    }

    // Report the final solution as the optimal. Candidates that tie the current best are accepted by DDS, so the final solution can differ from the
    // first parameter set found with the same objective.
    if (Fbest <= m_BestObjective) {
        m_BestObjective = Fbest;
        m_BestAlternative = bestParams;
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
void DiscreteDDSAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Iterations          : %d\n", m_MaxIter);
    fprintf(pFile, "Perturbation Value      : %lf\n", m_r_val);
    fprintf(pFile, "Initialization          : ");
    if (m_bWarmStart) { fprintf(pFile, "Warm Start\n"); }
    else if (m_UserSuppliedInit) { fprintf(pFile, "Initial Parameter Values\n"); }
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
void DiscreteDDSAlgorithm::WriteEndingMetrics(void) {

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
    fprintf(pFile, "Algorithm               : Discretely Dynamically-Dimensioned Search Algorithm (DDDS)\n");
    fprintf(pFile, "Desired Convergence Val : N/A\n");
    fprintf(pFile, "Actual Convergence Val  : N/A\n");
    fprintf(pFile, "Max Generations         : %d\n", m_MaxIter);
    fprintf(pFile, "Actual Generations      : %d\n", m_CurIter);
    fprintf(pFile, "Peterbation Value       : %lf\n", m_r_val);
    fprintf(pFile, "Num. Discrete Corrections : %d\n", m_nCorr);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);
    fprintf(pFile, "Algorithm successfully converged on a solution, however more runs may be needed\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
DiscreteDDS_Program()

Optimize the model using Discrete DDS.
**************************************************************************************************************************************************************/
void DiscreteDDS_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        // Create the Discrete DDS algorithm
        DiscreteDDSAlgorithm* DDDS = new DiscreteDDSAlgorithm();

        // Begin the primary worker operations. This transfers information to the secondary workers
        // and manages the solution process
        DDDS->Optimize();

    } else {
        // This is a secondary worker. Create the secondary worker and call the work routine
        // Create and setup the secondary worker, telling it to configure from MPI
        ModelWorker worker = ModelWorker(true);

        // Start the work in the processor
        worker.WorkMPI();

        // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
    }
} /* end DiscreteDDS_Program() */
