/**************************************************************************************************************************************************************

Simulated Annealing (SA) based on: Vanderbilt and Louie. 1984. "A Monte Carlo Simulated Annealing Approach to Optimization over Continuous Variables".
Journal of Computational Physics. vol. 56, pg. 259-271. The design space is first "melted" with a random walk to estimate the initial temperature. The system
is then cooled over a series of outer iterations. Within each outer iteration, a number of transitions are attempted and accepted according to the Metropolis
criterion. With the Vanderbilt-Louie transition method, the step distribution is adapted after each outer iteration using the covariance of the accepted
parameter sets. Transitions are generated in batches that are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm
base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "VandSA2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values. The step distribution
matrix is initialized to the identity matrix.
**************************************************************************************************************************************************************/
VandSA::VandSA() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    int numberOfParameters = m_pParamGroup->GetNumParams();
    IroncladString inputFilename = GetInFileName();

    // Read simulated annealing information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Simulated Annealing block within the file
        if (CheckToken(pFile, "BeginSimulatedAlg", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndSimulatedAlg", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginSimulatedAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndSimulatedAlg") == NULL) {

                if (strstr(line, "NumInitialTrials") != NULL) {
                    // Set the number of melting trials
                    sscanf(line, "%s %d", tmp, &m_NumMelts);
                }
                else if (strstr(line, "TemperatureScaleFactor") != NULL) {
                    // Set the temperature reduction factor
                    sscanf(line, "%s %lf", tmp, &m_TempFactor);
                }
                else if (strstr(line, "FinalTemperature") != NULL) {
                    // Set the final temperature
                    sscanf(line, "%s %lf", tmp, &m_FinalTemp);
                    m_bUserFinalTemp = true;
                }
                else if (strstr(line, "OuterIterations") != NULL) {
                    // Set the number of outer iterations
                    sscanf(line, "%s %d", tmp, &m_MaxOuter);
                }
                else if (strstr(line, "InnerIterations") != NULL) {
                    // Set the number of inner iterations
                    sscanf(line, "%s %d", tmp, &m_MaxInner);
                }
                else if (strstr(line, "ConvergenceVal") != NULL) {
                    // Set the convergence value
                    sscanf(line, "%s %lf", tmp, &m_StopVal);
                }
                else if (strstr(line, "TransitionMethod") != NULL) {
                    // Set the transition method
                    sscanf(line, "%s %s", tmp, tmp2);
                    MyStrLwr(tmp2);
                    if (strcmp(tmp2, "vanderbilt") == 0) { m_TransitionMethod = TRANS_VANDR; }
                    else if (strcmp(tmp2, "gauss") == 0) { m_TransitionMethod = TRANS_GAUSS; }
                }

                // Get the next line
                line = GetNxtDataLine(pFile, inputFilename);
            }
        }
        else {
            // The simulated annealing algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Initialize the matrices, with the step distribution matrix set to the identity matrix
    m_Q = std::vector<std::vector<double>>(numberOfParameters, std::vector<double>(numberOfParameters, 0.00));
    m_QT = m_Q;
    m_Cov = m_Q;
    m_Shape = m_Q;
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        m_Q[entryParameter][entryParameter] = 1.00;
    }

    // TODO: port support for the TelescopingStrategy option. Telescoping corrections are applied to each transition by the pre-refactor algorithm,
    // but telescoping has not yet been re-enabled within the Algorithm base class.

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the SA and it's member variables
**************************************************************************************************************************************************************/
void VandSA::Destroy(void) {

    m_Current.clear();
    m_Q.clear();
    m_QT.clear();
    m_Cov.clear();
    m_Shape.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void VandSA::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
CorrectParameters()

Applies any user-specified parameter corrections to a parameter set and returns the corrected values.
**************************************************************************************************************************************************************/
std::vector<double> VandSA::CorrectParameters(std::vector<double> parameters) {

    SetParameterGroup(parameters);
    PerformParameterCorrections();

    std::vector<double> corrected;
    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        corrected.push_back(m_pParamGroup->GetParamPtr(entryParameter)->GetEstimatedValueTransformed());
    }

    return corrected;
}

/**************************************************************************************************************************************************************
Initialize()

Evaluates the starting point along with one perturbation (+/-5% of the parameter range) for each additional MPI process, and returns the best objective.
Following the pre-refactor algorithm, each perturbation is applied to the previous perturbation.
**************************************************************************************************************************************************************/
double VandSA::Initialize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // The starting point is the first parameter set
    std::vector<std::vector<double>> samples;
    samples.push_back(CorrectParameters(m_Current));

    // Perturb the starting point for each of the remaining processes
    std::vector<double> perturbed = m_Current;
    for (int entrySample = 1; entrySample < numberOfMpiProcesses; entrySample++) {
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            double lwr = m_pParamGroup->GetParamPtr(entryParameter)->GetLowerBoundTransformed();
            double upr = m_pParamGroup->GetParamPtr(entryParameter)->GetUpperBoundTransformed();

            // Random value from -5% to +5% of the range
            double r = (UniformRandom() - 0.5) * 0.1 * (upr - lwr);
            r += perturbed[entryParameter];

            // Enforce bounds
            if (r > upr) r = upr;
            if (r < lwr) r = lwr;

            perturbed[entryParameter] = r;
        }

        perturbed = CorrectParameters(perturbed);
        samples.push_back(perturbed);
    }

    // Evaluate the starting points
    std::vector<double> objectives(samples.size(), INFINITY);
    ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);

    // The starting point sets the initial best, which is replaced by any improving perturbation
    double bestObjective = objectives[0];
    m_Current = samples[0];
    for (int entrySample = 1; entrySample < (int)samples.size(); entrySample++) {
        if (objectives[entrySample] < bestObjective) {
            bestObjective = objectives[entrySample];
            m_Current = samples[entrySample];
        }
    }

    m_BestObjective = bestObjective;
    m_BestAlternative = m_Current;

    return bestObjective;
}

/**************************************************************************************************************************************************************
GenerateRandomMove()

Displaces each parameter by a random amount, uniformly distributed within +/-10% of the parameter range. Moves that leave the bounds are replaced by moving
half the distance to the bound.
**************************************************************************************************************************************************************/
void VandSA::GenerateRandomMove(std::vector<double>& parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        double upr = pParam->GetUpperBoundTransformed();
        double lwr = pParam->GetLowerBoundTransformed();
        double curVal = parameters[entryParameter];
        double width = (upr - lwr);
        if ((width >= 1.00) && (width < 10.00) && (strcmp(pParam->GetType(), "real") != 0)) width = 10.00;

        // Generate a random displacement
        double r = (((2.00 * ((double)MyRand() / (double)MY_RAND_MAX)) - 1.00) * width) / 10.00;
        double adjVal = curVal + r;

        // If out of bounds, move half the distance
        if (adjVal > upr) { adjVal = (upr + curVal) / 2.00; m_NumUprViols++; }
        if (adjVal < lwr) { adjVal = (curVal + lwr) / 2.00; m_NumLwrViols++; }

        parameters[entryParameter] = adjVal;
    }
}

/**************************************************************************************************************************************************************
Melt()

'Melts' the design space to determine the initial temperature. A random walk is taken from the current best solution, and the initial temperature is
computed as the standard deviation of the objectives along the walk. The steps of the walk are evaluated together.

Returns the best objective function value.
**************************************************************************************************************************************************************/
double VandSA::Melt(double bestObjective) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Generate the random walk starting from the current best solution
    std::vector<std::vector<double>> samples;
    std::vector<double> walk = m_Current;
    for (int entrySample = 0; entrySample < m_NumMelts; entrySample++) {
        GenerateRandomMove(walk);
        walk = CorrectParameters(walk);
        samples.push_back(walk);
    }

    // Evaluate the steps
    std::vector<double> melts(samples.size(), INFINITY);
    ManageSingleObjectiveIterations(samples, 0, numberOfParameters, melts);
    m_MeltCount += m_NumMelts;

    // Update the best solution, if a new one is found
    for (int entrySample = 0; entrySample < m_NumMelts; entrySample++) {
        if (melts[entrySample] < bestObjective) {
            bestObjective = melts[entrySample];
            m_Current = samples[entrySample];
        }
    }

    // Compute the current convergence value
    double median = CalcMedian(melts.data(), m_NumMelts);
    m_CurStop = fabs((median - bestObjective) / median);

    /*---------------------------------------------
    Check the initial temperature against the
    following constraints:
       NEARLY_ZERO <= final temperature <=  1.00
    If these constraints are violated, suggest user
    adjust the 'OuterIterations' parameter.
    ---------------------------------------------*/
    m_CurTemp = m_InitTemp = CalcStdDev(melts.data(), m_NumMelts, CENTRAL_TEND_PCTILE);

    // If user supplied final temp, compute corresponding reduction factor
    if (m_bUserFinalTemp == true) {
        m_TempFactor = pow((m_FinalTemp / m_InitTemp), 1.00 / (double)m_MaxOuter);
    }

    if ((m_TempFactor < 1.00) && (m_TempFactor > 0.00)) {
        double maxTemp = 1.00 / pow(m_TempFactor, (double)m_MaxOuter);
        double minTemp = NEARLY_ZERO * maxTemp;

        if (m_InitTemp < minTemp) {
            LogError(ERR_SA_TEMP, "Final temperature nearly zero, consider reducing OuterIterations");
        }
        if (m_InitTemp > maxTemp) {
            LogError(ERR_SA_TEMP, "Final temperature very high, consider increasing OuterIterations");
        }
    }
    else {
        // User-supplied reduction rate is not valid, compute internally
        LogError(ERR_BAD_ARGS, "Invalid temperature reduction rate; using internally calculated value");
        m_TempFactor = pow((1.00 / m_InitTemp), 1.00 / (double)m_MaxOuter);
    }

    // Track the best solution
    if (bestObjective < m_BestObjective) {
        m_BestObjective = bestObjective;
        m_BestAlternative = m_Current;
    }

    return bestObjective;
}

/**************************************************************************************************************************************************************
VanderbiltTransition()

Generates a move from the current parameter set using equation (5) of the paper: dx = Q*u, where u is uniformly distributed from -sqrt(3) to +sqrt(3).
Moves that leave the bounds are replaced by moving half the distance to the bound.
**************************************************************************************************************************************************************/
std::vector<double> VandSA::VanderbiltTransition(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    std::vector<double> proposal = m_Current;

    // Randomly initialize the 'u' vector
    std::vector<double> u;
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        u.push_back((((double)MyRand() / (double)MY_RAND_MAX) * 2.00 * sqrt(3.00)) - sqrt(3.00));
    }

    // Compute delta x
    std::vector<double> dx(numberOfParameters, 0.00);
    std::vector<double*> Q;
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) { Q.push_back(m_Q[entryParameter].data()); }
    VectMult(Q.data(), u.data(), dx.data(), numberOfParameters, numberOfParameters);

    // Make a move
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        double upr = pParam->GetUpperBoundTransformed();
        double lwr = pParam->GetLowerBoundTransformed();
        double val = m_Current[entryParameter];
        double adj = val + dx[entryParameter];

        // If out of bounds, move half the distance
        if (adj > upr) { adj = (upr + val) / 2.00; m_NumUprViols++; }
        if (adj < lwr) { adj = (val + lwr) / 2.00; m_NumLwrViols++; }

        proposal[entryParameter] = adj;
    }

    return CorrectParameters(proposal);
}

/**************************************************************************************************************************************************************
GaussTransition()

Generates a move from the current parameter set. Each parameter is perturbed using a normal distribution centered on the current value, with a standard
deviation estimated from the current objective function value (limited to 68% of the parameter range).
**************************************************************************************************************************************************************/
std::vector<double> VandSA::GaussTransition(double currentObjective) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    std::vector<double> proposal = m_Current;
    double sd = sqrt(MyMax(NEARLY_ZERO, fabs(currentObjective)) / (double)numberOfParameters);

    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        double upr = pParam->GetUpperBoundTransformed();
        double lwr = pParam->GetLowerBoundTransformed();
        double curVal = m_Current[entryParameter];

        // 1 std. dev. should not cover more than 68% of the range
        double range = upr - lwr;
        double sdmax = range * 0.68;
        double sdi = (sd > sdmax) ? sdmax : sd;

        // Epsilon perturbation using normal distribution
        double val = MyGaussRand(curVal, sdi);

        // Enforce parameter limits
        double r = 2 * ((double)MyRand() / (double)MY_RAND_MAX) - 1.00; // -1 to 1
        if ((val > upr) && (r >= 0.00)) val = curVal + (upr - curVal) * r;
        if ((val > upr) && (r < 0.00)) val = curVal + (curVal - lwr) * r;
        if ((val < lwr) && (r >= 0.00)) val = curVal - (curVal - lwr) * r;
        if ((val < lwr) && (r < 0.00)) val = curVal - (upr - curVal) * r;

        proposal[entryParameter] = val;
    }

    return CorrectParameters(proposal);
}

/**************************************************************************************************************************************************************
AcceptTransition()

Tests a solved transition against the acceptance criteria: the move is accepted if it does not increase the objective, or if a uniform random number is less
than the acceptance probability exp(-dE/T). Accepted moves become the current parameter set. Returns the objective function value at the current location.
**************************************************************************************************************************************************************/
double VandSA::AcceptTransition(std::vector<double> proposal, double proposalObjective, double currentObjective) {

    if (proposalObjective <= currentObjective) {
        // Accept the downhill move
        m_NumDownhill++;
        m_Current = proposal;
        return proposalObjective;
    }

    // Test against acceptance probability
    double increase = proposalObjective - currentObjective;
    double r = (double)MyRand() / (double)MY_RAND_MAX;
    double prob = exp(-increase / m_CurTemp);

    // Probability metrics
    m_TotProb += prob;
    m_NumProbTests++;

    if (prob >= r) {
        // Accept the uphill move
        m_NumUphill++;
        m_Current = proposal;
        return proposalObjective;
    }

    // Reject the move
    m_NumAborts++;
    return currentObjective;
}

/**************************************************************************************************************************************************************
UpdateStepMatrix()

Updates the step distribution matrix (Q) from the parameter sets visited during an equilibration. The average parameter values (eqn 10 of the paper) and
shape matrix (eqn 11) are used to estimate the covariance of the next iteration (eqn 13), which is factored using Cholesky decomposition (eqn 7).
**************************************************************************************************************************************************************/
void VandSA::UpdateStepMatrix(std::vector<std::vector<double>> transitions) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Average parameter transition (eqn 10 of the paper)
    std::vector<double> A(numberOfParameters, 0.00);
    for (int i = 0; i < numberOfParameters; i++) {
        for (int m = 0; m < m_MaxInner; m++) {
            A[i] += transitions[m][i];
        }
        A[i] /= (double)m_MaxInner;
    }

    // Shape estimate (eqn 11 of the paper)
    for (int i = 0; i < numberOfParameters; i++) {
        for (int j = 0; j < numberOfParameters; j++) {
            m_Shape[i][j] = 0.00;
            for (int m = 0; m < m_MaxInner; m++) {
                m_Shape[i][j] += (transitions[m][i] - A[i]) * (transitions[m][j] - A[j]);
            }
            m_Shape[i][j] /= (double)m_MaxInner;
        }
    }

    // Estimate the covariance of the next iteration (eqn 13 of the paper)
    for (int i = 0; i < numberOfParameters; i++) {
        for (int j = 0; j < numberOfParameters; j++) {
            m_Cov[i][j] = (3.00 * m_Shape[i][j]) / (0.11 * (double)m_MaxInner);
        }
    }

    // Compute Q of the next iteration, using Cholesky decomposition (eqn 7 of the paper)
    std::vector<double*> cov, Q, QT;
    for (int i = 0; i < numberOfParameters; i++) {
        cov.push_back(m_Cov[i].data());
        Q.push_back(m_Q[i].data());
        QT.push_back(m_QT[i].data());
    }
    CholeskyDecomp(cov.data(), Q.data(), QT.data(), numberOfParameters);
}

/**************************************************************************************************************************************************************
Equilibrate()

Allows the system to come to equilibrium at the current temperature. This is accomplished by making a number of moves equal to m_MaxInner, starting from the
best solution. Moves are generated in batches from the current location, solved together, and then tested for acceptance in order. The parameter sets
visited are used to update the step distribution matrix. Returns the objective function value of the best move.
**************************************************************************************************************************************************************/
double VandSA::Equilibrate(double bestObjective) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    double currentObjective = bestObjective;
    std::vector<double> finner;                                             // Current objective after each transition (used in convergence test)
    std::vector<std::vector<double>> transitions;                           // Current parameter set after each transition

    // Initialize probability metrics for this set of transitions
    m_NumProbTests = 0;
    m_TotProb = 0.00;

    while ((int)finner.size() < m_MaxInner) {
        // Generate a batch of transitions from the current location
        int batchSize = std::min(m_BatchSize, m_MaxInner - (int)finner.size());
        std::vector<std::vector<double>> samples;
        for (int entrySample = 0; entrySample < batchSize; entrySample++) {
            if (m_TransitionMethod == TRANS_VANDR) {
                samples.push_back(VanderbiltTransition());
            }
            else {
                samples.push_back(GaussTransition(currentObjective));
            }
            m_TransCount++;
        }

        // Evaluate the transitions
        std::vector<double> objectives(samples.size(), INFINITY);
        ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);

        // Accept or reject each transition in order
        for (int entrySample = 0; entrySample < batchSize; entrySample++) {
            currentObjective = AcceptTransition(samples[entrySample], objectives[entrySample], currentObjective);
            finner.push_back(currentObjective);
            transitions.push_back(m_Current);

            // Update the best solution
            if (currentObjective < bestObjective) {
                bestObjective = currentObjective;
                m_BestObjective = currentObjective;
                m_BestAlternative = m_Current;
            }
        }
    }

    // Compute the current convergence value (eqn 18 of the paper)
    double median = CalcMedian(finner.data(), m_MaxInner);
    m_CurStop = fabs((median - bestObjective) / median);

    // Update the step distribution matrix
    UpdateStepMatrix(transitions);

    // Compute avg. probability of acceptance for the equilibration
    if (m_NumProbTests > 0) {
        m_CurProb = m_TotProb / (double)m_NumProbTests;
        if (m_InitProb < 0.00) { m_InitProb = m_CurProb; }
    }

    // Subsequent transitions start from the best solution
    m_Current = m_BestAlternative;
    m_EquilCount++;

    return bestObjective;
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using Vanderbilt-Louie SA.
**************************************************************************************************************************************************************/
void VandSA::Optimize(void) {

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

    // Transitions are evaluated in batches, one per available worker
    m_BatchSize = numberOfMpiProcesses - 1;
    if (m_bSolveOnPrimary) { m_BatchSize++; }
    if (m_BatchSize < 1) { m_BatchSize = 1; }

    if (m_bWarmStart && (int)m_BestAlternative.size() >= numberOfParameters) {
        // Start from the best solution of the previous run
        m_Current = std::vector<double>(m_BestAlternative.begin(), m_BestAlternative.begin() + numberOfParameters);
    }
    else {
        // Start from the initial parameter values
        m_Current.clear();
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            m_Current.push_back(m_pParamGroup->GetParamPtr(entryParameter)->GetInitialValueTransformed());
        }
    }

    if (!m_bWarmStart) {
        // Start a clean analysis
        // Write the information to the primary log file
        WriteSetup2("VandSA", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();
    }

    // Evaluate the starting point and melt the design space
    double bestObjective = Initialize();
    bestObjective = Melt(bestObjective);

    // Write the initial result to the log file
    m_NumOuter = 0;
    SetParameterGroup(m_BestAlternative);
    WriteRecord2(m_NumOuter, m_BestObjective, m_CurStop, m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Enter the cooling loop
    for (int entryOuter = 0; entryOuter < m_MaxOuter; entryOuter++) {

        bestObjective = Equilibrate(bestObjective);
        m_CurTemp *= m_TempFactor;
        m_NumOuter = entryOuter + 1;

        // Write the iteration result to the log file
        SetParameterGroup(m_BestAlternative);
        WriteRecord2(m_NumOuter, m_BestObjective, m_CurStop, m_pObsGroup, m_pObjFunc, m_pParamGroup);

        // Check for convergence
        if (m_CurStop <= m_StopVal) { break; }
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
void VandSA::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
    fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
    fprintf(pFile, "Melting Trials          : %d\n", m_NumMelts);
    fprintf(pFile, "Final Temperature       : ");
    if (m_bUserFinalTemp) { fprintf(pFile, "%E\n", m_FinalTemp); }
    else { fprintf(pFile, "Temperature Reduction Factor of %.2lf\n", m_TempFactor); }
    fprintf(pFile, "Transition Method       : ");
    if (m_TransitionMethod == TRANS_VANDR) { fprintf(pFile, "Vanderbilt-Louie\n"); }
    else { fprintf(pFile, "Gaussian\n"); }

    // Write the parameter header to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Convergence\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void VandSA::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Simulated Annealing (Vanderbilt-Louie Implementation)\n");
    fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
    fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
    fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
    fprintf(pFile, "Actual Outer Iterations : %d\n", m_NumOuter);
    fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
    fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor * 100.0);
    fprintf(pFile, "Initial Temperature     : %E\n", m_InitTemp);
    fprintf(pFile, "Final Temperature       : %E\n", m_CurTemp);
    fprintf(pFile, "Initial Pr[Acc]         : %.2lf%%\n", m_InitProb * 100.0);
    fprintf(pFile, "Actual Final Pr[Acc]    : %.2lf%%\n", m_CurProb * 100.0);
    fprintf(pFile, "Expected Final Pr[Acc]  : 50.00%%\n");
    fprintf(pFile, "Melting Evals           : %d\n", m_MeltCount);
    if (m_TransitionMethod == TRANS_VANDR) { fprintf(pFile, "Transition Method       : Vanderbilt-Louie\n"); }
    else { fprintf(pFile, "Transition Method       : Gaussian\n"); }
    fprintf(pFile, "Transition Evals        : %d\n", m_TransCount);
    fprintf(pFile, "Equilibration Evals     : %d\n", m_EquilCount);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);
    fprintf(pFile, "Rejected Transitions    : %d\n", m_NumAborts);
    fprintf(pFile, "Uphill Transitions      : %d\n", m_NumUphill);
    fprintf(pFile, "Downhill Transitions    : %d\n", m_NumDownhill);
    fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
    fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

    if (m_CurStop <= m_StopVal) {
        fprintf(pFile, "Algorithm successfully converged on a solution\n");
    }
    else {
        fprintf(pFile, "Algorithm failed to converge on a solution, more outer iterations may be needed\n");
    }

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
VSA_Program()

Optimize the model using Vanderbilt-Louie SA.
**************************************************************************************************************************************************************/
void VSA_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the Vanderbilt-Louie simulated annealing algorithm
       VandSA *VSA = new VandSA();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       VSA->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
} /* end VSA_Program() */
