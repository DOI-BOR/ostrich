/**************************************************************************************************************************************************************

Simulated Annealing (SA) for combinatorial and integer (discrete) parameters. The design space is first "melted" with a random walk to estimate the initial
temperature. The system is then cooled over a series of outer iterations. Within each outer iteration, a number of transitions (moves of -1, 0 or +1 for
each parameter) are attempted and are accepted if they improve the objective or, otherwise, with a probability that decreases as the temperature drops
(Metropolis criterion). Transitions are generated in batches that are solved using the shared, asynchronous MPI work-farming functionality provided by the
Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "ComboSA2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. All parameters are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
ComboSA::ComboSA() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
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

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the SA and it's member variables
**************************************************************************************************************************************************************/
void ComboSA::Destroy(void) {

    m_Current.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Sets a parameter set into the parameter group so that it is reflected in the log/model state.
**************************************************************************************************************************************************************/
void ComboSA::SetParameterGroup(std::vector<double> parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
GenerateRandomMove()

Displaces each parameter by -1, 0 or +1 (chosen at random), wrapping around the parameter bounds.
**************************************************************************************************************************************************************/
void ComboSA::GenerateRandomMove(std::vector<double>& parameters) {

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        int upr = (int)(pParam->GetUpperBoundTransformed());
        int lwr = (int)(pParam->GetLowerBoundTransformed());
        int val = (int)(parameters[entryParameter] - lwr);
        int range = (upr - lwr);

        // Generate a random displacement
        int adj = (MyRand() % 3) - 1; // -1, 0, +1

        // Perturb parameter value
        val += adj;
        if (val < 0) val += (range + 1);
        else         val %= (range + 1);

        // If out of bounds (shouldn't happen, but check anyway)
        if (val > range) { val = range; m_NumUprViols++; }
        if (val < 0) { val = 0; m_NumLwrViols++; }

        parameters[entryParameter] = (double)(val + lwr);
    }
}

/**************************************************************************************************************************************************************
Melt()

'Melts' the design space to determine the initial temperature. A random walk is taken from the starting point, and the energy change (dE) between consecutive
steps of the walk is recorded. The initial temperature is computed so that a statistically large energy increase will be accepted with ~100% probability.
That is:
   T0 = -dEmax/ln(0.99) ~ 100*dEmax
Where dEmax is computed as 3 std. deviations from the mean dE estimated from the melting trials. The steps of the walk are evaluated together.

Returns the best objective function value.
**************************************************************************************************************************************************************/
double ComboSA::Melt(double bestObjective) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Generate the random walk starting from the current solution
    std::vector<std::vector<double>> samples;
    std::vector<double> walk = m_Current;
    for (int entrySample = 0; entrySample < m_NumMelts; entrySample++) {
        GenerateRandomMove(walk);
        samples.push_back(walk);
    }

    // Evaluate the steps
    std::vector<double> objectives(samples.size(), INFINITY);
    ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);

    // Process the steps in the order of the walk
    std::vector<double> melts;                                              // Objectives of the melting trials
    std::vector<double> pdE;                                                // Magnitude of the energy changes
    double previousObjective = bestObjective;
    double dEavg = 0.00;

    for (int entrySample = 0; entrySample < m_NumMelts; entrySample++) {
        double objective = objectives[entrySample];
        double dE = objective - previousObjective;
        previousObjective = objective;
        m_MeltCount++;

        // Accumulate values
        melts.push_back(objective);
        pdE.push_back(fabs(dE));
        dEavg += fabs(dE);

        // Decrease in energy, check for a new best solution
        if ((dE < 0.00) && (objective < bestObjective)) {
            bestObjective = objective;
            m_Current = samples[entrySample];
        }
    }

    dEavg /= (double)m_NumMelts;

    // Compute the current convergence value
    double Emed = CalcMedian(melts.data(), m_NumMelts);
    m_CurStop = fabs((Emed - bestObjective) / Emed);

    // Assign the initial temperature
    double dEmax = dEavg + 3.00 * CalcStdDev(pdE.data(), m_NumMelts, CENTRAL_TEND_PCTILE);
    m_CurTemp = m_InitTemp = 100.00 * dEmax;
    m_dEavg = dEavg; // store average energy change as a metric

    // If user-supplied reduction rate is not valid, compute internally
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
Transition()

Generates a move from the current parameter set. Each parameter is perturbed to a neighboring state (-1, 0 or +1), wrapping around the parameter bounds.
**************************************************************************************************************************************************************/
std::vector<double> ComboSA::Transition(void) {

    std::vector<double> proposal = m_Current;

    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        ParameterABC* pParam = m_pParamGroup->GetParamPtr(entryParameter);
        int upr = (int)(pParam->GetUpperBoundTransformed());
        int lwr = (int)(pParam->GetLowerBoundTransformed());
        int range = upr - lwr;
        int val = (int)(m_Current[entryParameter]) - lwr;

        // Determine random direction
        int dir = (MyRand() % 3) - 1;

        // Perturb parameter value
        val += dir;
        if (val < 0) val += (range + 1);
        else         val %= (range + 1);

        // Modulus arithmetic should prevent bounds violations
        if (val > range) { val = range; m_NumUprViols++; }
        if (val < 0) { val = 0; m_NumLwrViols++; }

        proposal[entryParameter] = (double)(val + lwr);
    }

    return proposal;
}

/**************************************************************************************************************************************************************
AcceptTransition()

Tests a solved transition against the acceptance criteria: the move is accepted if it decreases the objective, or if a uniform random number is less than
the acceptance probability exp(-dE/T). Accepted moves become the current parameter set. Returns the objective function value at the current location.
**************************************************************************************************************************************************************/
double ComboSA::AcceptTransition(std::vector<double> proposal, double proposalObjective, double currentObjective) {

    if (proposalObjective < currentObjective) {
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
double ComboSA::Equilibrate(double bestObjective) {

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
            samples.push_back(Transition());
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

Minimize the objective function using discrete SA.
**************************************************************************************************************************************************************/
void ComboSA::Optimize(void) {

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
        WriteSetup2("ComboSA", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();
    }

    // Evaluate the starting point
    std::vector<std::vector<double>> samples = { m_Current };
    std::vector<double> objectives(1, INFINITY);
    ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);
    m_MeltCount++;

    double bestObjective = objectives[0];
    m_BestObjective = bestObjective;
    m_BestAlternative = m_Current;

    // Melt the design space
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
void ComboSA::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
    fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
    fprintf(pFile, "Melting Trials          : %d\n", m_NumMelts);
    fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor * 100.0);

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
void ComboSA::WriteEndingMetrics(void) {

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
    fprintf(pFile, "Algorithm               : Simulated Annealing for Discrete Parameters\n");
    fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
    fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
    fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
    fprintf(pFile, "Actual Outer Iterations : %d\n", m_NumOuter);
    fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
    fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor * 100.0);
    fprintf(pFile, "Initial Temperature     : %.2lf\n", m_InitTemp);
    fprintf(pFile, "Avg. Energy Change      : %.2lf\n", m_dEavg);
    fprintf(pFile, "Std. Dev. Energy Change : %.2lf\n", sd);
    fprintf(pFile, "Final Temperature       : %.2lf\n", m_CurTemp);
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
CSA_Program()

Optimize the model using discrete (combinatorial) SA.
**************************************************************************************************************************************************************/
void CSA_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the discrete simulated annealing algorithm
       ComboSA *CSA = new ComboSA();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       CSA->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
} /* end CSA_Program() */
