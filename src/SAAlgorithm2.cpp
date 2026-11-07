/**************************************************************************************************************************************************************

Simulated Annealing (SA) for continuous parameters. The design space is first "melted" with a random walk to estimate the initial temperature. The system is
then cooled over a series of outer iterations. Within each outer iteration, a number of transitions (moves to neighboring solutions) are attempted and are
accepted if they improve the objective or, otherwise, with a probability that decreases as the temperature drops (Metropolis criterion). Transitions are
generated in batches that are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "SAAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
SAAlgorithm::SAAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
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
                    // Set the final temperature, or the method used to compute it
                    sscanf(line, "%s %s", tmp, tmp2);
                    MyStrLwr(tmp2);
                    if ((strcmp(tmp2, "computed-vanderbilt") == 0) || (strcmp(tmp2, "computed") == 0)) { m_FinalTempMethod = TMETHD_VNDR; }
                    else if (strcmp(tmp2, "computed-ben-ameur") == 0) { m_FinalTempMethod = TMETHD_BAMR; }
                    else {
                        m_FinalTempMethod = TMETHD_USER;
                        m_FinalTemp = atof(tmp2);
                    }
                }
                else if (strstr(line, "TransitionMethod") != NULL) {
                    // Set the transition method
                    sscanf(line, "%s %s", tmp, tmp2);
                    MyStrLwr(tmp2);
                    if (strcmp(tmp2, "uniform") == 0) { m_TransitionMethod = TRANS_UNFRM; }
                    else if (strcmp(tmp2, "gauss") == 0) { m_TransitionMethod = TRANS_GAUSS; }
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

    // TODO: port support for the TelescopingStrategy option. Telescoping corrections are applied to each transition by the pre-refactor algorithm,
    // but telescoping has not yet been re-enabled within the Algorithm base class.

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the SA and it's member variables
**************************************************************************************************************************************************************/
void SAAlgorithm::Destroy(void) {

    m_Current.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void SAAlgorithm::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
CorrectParameters()

Applies any user-specified parameter corrections to a parameter set and returns the corrected values.
**************************************************************************************************************************************************************/
std::vector<double> SAAlgorithm::CorrectParameters(std::vector<double> parameters) {

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
double SAAlgorithm::Initialize(void) {

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

Displaces each parameter by a random amount, uniformly distributed within a 20% window (+/-10% of the range) around the current value. The window is shifted
away from the bounds as needed.
**************************************************************************************************************************************************************/
void SAAlgorithm::GenerateRandomMove(std::vector<double>& parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        double curVal = parameters[entryParameter];
        double upr = pParam->GetUpperBoundTransformed();
        double lwr = pParam->GetLowerBoundTransformed();
        double range = upr - lwr;

        if ((range >= 1.00) && (range < 10.00) && (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;

        // +/-10% of the range = range/5
        double range10pct = range / 5;

        // Uniform random number between 0 and 1
        double r = ((double)MyRand() / (double)MY_RAND_MAX);

        // Minimum value after perturbation --- logic handles edge cases
        double vmin, vmax;
        double vavg = (upr + lwr) * 0.5;
        if (curVal < vavg) {
            vmin = curVal - 0.5 * range10pct;
            if (vmin < lwr) vmin = lwr;
            vmax = vmin + range10pct;
        }
        else {
            vmax = curVal + 0.5 * range10pct;
            if (vmax > upr) vmax = upr;
            vmin = vmax - range10pct;
        }

        // Adjusted value = (10% range * random) + minimum
        double adjVal = (range10pct * r) + vmin;

        // Enforce limits by moving half the distance --- this shouldn't happen given logic above!
        if (adjVal > upr) { adjVal = (upr + curVal) / 2.00; m_NumUprViols++; }
        if (adjVal < lwr) { adjVal = (curVal + lwr) / 2.00; m_NumLwrViols++; }

        parameters[entryParameter] = adjVal;
    }
}

/**************************************************************************************************************************************************************
Melt()

'Melts' the design space to determine the initial temperature. A random walk is taken from the current best solution, and the energy change (dE) between
consecutive steps of the walk is recorded. The initial temperature is computed so that a statistically large energy increase will be accepted with ~100%
probability. That is:
   T0 = -dEmax/ln(0.99) ~ 100*dEmax
Where dEmax is computed as 2 std. deviations from the mean dE estimated from the melting trials. The steps of the walk are evaluated together, and any step
that overflows is discarded and replaced by continuing the walk.

Returns the best objective function value.
**************************************************************************************************************************************************************/
double SAAlgorithm::Melt(double bestObjective) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    std::vector<double> melts;                                              // Objectives of the melting trials
    std::vector<double> pdE;                                                // Magnitude of the energy changes
    double previousObjective = bestObjective;
    double dEavg = 0.00;

    // Random walk starting from the current best solution
    std::vector<double> walk = m_Current;

    while ((int)melts.size() < m_NumMelts) {
        // Generate the next set of steps in the walk
        std::vector<std::vector<double>> samples;
        for (int entrySample = (int)melts.size(); entrySample < m_NumMelts; entrySample++) {
            GenerateRandomMove(walk);
            walk = CorrectParameters(walk);
            samples.push_back(walk);
        }

        // Evaluate the steps
        std::vector<double> objectives(samples.size(), INFINITY);
        ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);

        // Process the steps in the order of the walk
        for (int entrySample = 0; entrySample < (int)samples.size(); entrySample++) {
            double objective = objectives[entrySample];
            m_MeltCount++;

            // Discard results that overflow
            if (CheckOverflow(objective * objective) == true) {
                continue;
            }

            double dE = objective - previousObjective;
            previousObjective = objective;

            melts.push_back(objective);
            pdE.push_back(fabs(dE));
            dEavg += fabs(dE);

            // Update the best solution, if a new one is found
            if (objective < bestObjective) {
                bestObjective = objective;
                m_Current = samples[entrySample];
            }
        }
    }

    // Average or median energy change, use whichever is smaller to avoid influence of outliers
    dEavg /= (double)m_NumMelts;
    double dEmed = CalcMedian(pdE.data(), m_NumMelts);
    if (dEmed < dEavg) dEavg = dEmed;

    // Compute the current convergence value
    double Emed = CalcMedian(melts.data(), m_NumMelts);
    m_CurStop = fabs((Emed - bestObjective) / Emed);

    // Assign the initial temperature --- shoot for 95th percentile of energy change
    double dEmax = dEavg + 2.00 * CalcStdDev(pdE.data(), m_NumMelts, CENTRAL_TEND_PCTILE);
    m_CurTemp = m_InitTemp = 100.00 * dEmax;
    m_dEavg = dEavg; // store average energy change as a metric

    if (m_FinalTempMethod == TMETHD_USER) {
        // User supplied final temp, compute corresponding reduction factor
        m_TempFactor = pow((m_FinalTemp / m_InitTemp), 1.00 / (double)m_MaxOuter);
    }
    else if (m_FinalTempMethod == TMETHD_BAMR) {
        /* -----------------------------------------------------------------
        Use algorithm described by Walid Ben-Ameur in "Computing the Initial
        Temperature of Simulated Annealing", 2004, Computational Optimization
        and Applications, vol. 29, pg. 369-385.
        ----------------------------------------------------------------- */
        // Step 1a -- number of transitions is one less than melting trials
        int S = (m_NumMelts - 1);
        // Step 1b -- list of positive transitions from melting trials
        double* pS = melts.data();
        // Step 1c -- Set T1 equal to previously computed init. temperature
        double T1 = m_InitTemp;

        // Step 2 -- iteratively compute Pn (the acceptance prob.) and Tn. Repeated for the initial (P0 = 0.99) and final (P0 = 0.01) temperatures
        double P0Values[2] = { 0.99, 0.01 };
        double* pTValues[2] = { &m_InitTemp, &m_FinalTemp };

        for (int entryPass = 0; entryPass < 2; entryPass++) {
            double P0 = P0Values[entryPass];
            double dpLast = 2.00;
            double p = 1.00;
            double Tn = T1;

            for (int t = 0; t < 10000; t++) {
                double numer = 0.00;
                double denom = 0.00;
                for (int i = 0; i < S; i++) {
                    double Emax, Emin;
                    if (pS[i] > pS[i + 1]) { Emax = pS[i]; Emin = pS[i + 1]; }
                    else { Emin = pS[i]; Emax = pS[i + 1]; }
                    numer += exp(-Emax / Tn);
                    denom += exp(-Emin / Tn);
                }
                double Pn = numer / denom;
                double dP = fabs(Pn - P0);
                if (dP <= 0.001) break;

                // Divergence check
                if (dP >= dpLast) {
                    p *= 2.00;
                }
                else {
                    Tn = Tn * pow(log(Pn) / log(P0), (1.00 / p));
                    dpLast = dP;
                }
            }
            *pTValues[entryPass] = Tn;
        }

        m_TempFactor = pow((m_FinalTemp / m_InitTemp), 1.00 / (double)m_MaxOuter);
        m_CurTemp = m_InitTemp;
    }
    else if (m_FinalTempMethod == TMETHD_VNDR) {
        // Pre-computed initial and final temperatures
        m_InitTemp = -m_dEavg / log(0.99);
        m_FinalTemp = -m_dEavg / log(0.01);
        m_TempFactor = pow((m_FinalTemp / m_InitTemp), 1.00 / (double)m_MaxOuter);
        m_CurTemp = m_InitTemp;
    }

    // If reduction rate is not valid, compute internally
    if ((m_TempFactor >= 1.00) || (m_TempFactor <= 0.00)) {
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
UniformTransition()

Generates a move from the current parameter set. Each parameter is perturbed within a 20% window (+/-10%) of either the parameter range or the current
value (chosen at random for each move).
**************************************************************************************************************************************************************/
std::vector<double> SAAlgorithm::UniformTransition(void) {

    std::vector<double> proposal = m_Current;

    // Make a move --- randomly perturb parameters relative to range or relative to current value
    unsigned int rangeType = MyRand() % 2;
    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        double curval = m_Current[entryParameter];
        double upr = pParam->GetUpperBoundTransformed();
        double lwr = pParam->GetLowerBoundTransformed();
        double range;

        if (rangeType == 0) {
            // Absolute range
            range = upr - lwr;
        }
        else {
            // Relative range
            range = curval;
            if (range <= NEARLY_ZERO) {
                range = upr - lwr;
            }
        }

        if ((range >= 1.00) && (range < 10.00) && (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;

        // +/-10% of the range = range/5
        double range10pct = range / 5;

        // Uniform random number between 0 and 1
        double r = ((double)MyRand() / (double)MY_RAND_MAX);

        // Minimum value after perturbation --- logic handles edge cases
        double vmin, vmax;
        double vavg = (upr + lwr) * 0.5;
        if (curval < vavg) {
            vmin = curval - 0.5 * range10pct;
            if (vmin < lwr) vmin = lwr;
            vmax = vmin + range10pct;
        }
        else {
            vmax = curval + 0.5 * range10pct;
            if (vmax > upr) vmax = upr;
            vmin = vmax - range10pct;
        }

        // Adjusted value = (10% range * random) + minimum
        double val = (range10pct * r) + vmin;

        // Enforce limits by moving half the distance
        if (val > upr) { val = (upr + curval) / 2.00; m_NumUprViols++; }
        if (val < lwr) { val = (curval + lwr) / 2.00; m_NumLwrViols++; }

        proposal[entryParameter] = val;
    }

    return CorrectParameters(proposal);
}

/**************************************************************************************************************************************************************
GaussTransition()

Generates a move from the current parameter set. Each parameter is perturbed using a normal distribution centered on the current value, with a standard
deviation estimated from the current objective function value (limited to 68% of the parameter range).
**************************************************************************************************************************************************************/
std::vector<double> SAAlgorithm::GaussTransition(double currentObjective) {

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
double SAAlgorithm::AcceptTransition(std::vector<double> proposal, double proposalObjective, double currentObjective) {

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
Equilibrate()

Allows the system to come to equilibrium at the current temperature. This is accomplished by making a number of moves equal to m_MaxInner, starting from the
best solution. Moves are generated in batches from the current location, solved together, and then tested for acceptance in order. Returns the objective
function value of the best move.
**************************************************************************************************************************************************************/
double SAAlgorithm::Equilibrate(double bestObjective) {

    int numberOfParameters = m_pParamGroup->GetNumParams();
    double currentObjective = bestObjective;
    std::vector<double> finner;                                             // Current objective after each transition (used in convergence test)

    // Initialize probability metrics for this set of transitions
    m_NumProbTests = 0;
    m_TotProb = 0.00;

    while ((int)finner.size() < m_MaxInner) {
        // Generate a batch of transitions from the current location
        int batchSize = std::min(m_BatchSize, m_MaxInner - (int)finner.size());
        std::vector<std::vector<double>> samples;
        for (int entrySample = 0; entrySample < batchSize; entrySample++) {
            if (m_TransitionMethod == TRANS_UNFRM) {
                samples.push_back(UniformTransition());
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

            // Update the best solution
            if (currentObjective < bestObjective) {
                bestObjective = currentObjective;
                m_BestObjective = currentObjective;
                m_BestAlternative = m_Current;
            }
        }
    }

    // Compute the current convergence value
    double median = CalcMedian(finner.data(), m_MaxInner);
    m_CurStop = fabs((median - bestObjective) / median);

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

Minimize the objective function using SA.
**************************************************************************************************************************************************************/
void SAAlgorithm::Optimize(void) {

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
        WriteSetup2("SAAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
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
void SAAlgorithm::WriteStartingMetrics(void) {

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
    if (m_FinalTempMethod == TMETHD_USER) { fprintf(pFile, "%E\n", m_FinalTemp); }
    else if (m_FinalTempMethod == TMETHD_VNDR) { fprintf(pFile, "Computed (Vanderbilt)\n"); }
    else if (m_FinalTempMethod == TMETHD_BAMR) { fprintf(pFile, "Computed (Ben-Ameur)\n"); }
    else { fprintf(pFile, "Temperature Reduction Factor of %.2lf\n", m_TempFactor); }
    fprintf(pFile, "Transition Method       : ");
    if (m_TransitionMethod == TRANS_UNFRM) { fprintf(pFile, "Uniform\n"); }
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
void SAAlgorithm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    double sd = ((m_InitTemp / 100.00) - m_dEavg) / 3.00;
    double p = exp(-m_dEavg / m_CurTemp);
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Simulated Annealing for Continuous Parameters\n");
    fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
    fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
    fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
    fprintf(pFile, "Actual Outer Iterations : %d\n", m_NumOuter);
    fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
    fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor * 100.0);
    fprintf(pFile, "Initial Temperature     : %E\n", m_InitTemp);
    fprintf(pFile, "Avg. Energy Change      : %E\n", m_dEavg);
    fprintf(pFile, "Std. Dev. Energy Change : %E\n", sd);
    fprintf(pFile, "Final Temperature       : %E\n", m_CurTemp);
    fprintf(pFile, "Initial Pr[Acc]         : %.2lf%%\n", m_InitProb * 100.0);
    fprintf(pFile, "Actual Final Pr[Acc]    : %.2lf%%\n", m_CurProb * 100.0);
    fprintf(pFile, "Expected Final Pr[Acc]  : %.2lf%%\n", p * 100.0);
    fprintf(pFile, "Melting Evals           : %d\n", m_MeltCount);
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
SA_Program()

Optimize the model using SA.
**************************************************************************************************************************************************************/
void SA_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the simulated annealing algorithm
       SAAlgorithm *SA = new SAAlgorithm();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       SA->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
} /* end SA_Program() */
