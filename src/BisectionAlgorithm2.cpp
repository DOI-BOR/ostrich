/**************************************************************************************************************************************************************

The Bisection algorithm optimizes one parameter at a time via simple bisection (a coordinate-descent style search). For each parameter, the algorithm 
subdivides the parameter's bounded range into quadrants, evaluates the objective function at the resulting quadrant boundaries, and repeatedly discards 
the quadrant(s) least likely to contain the optimum. Parameters are processed one at a time (holding all other parameters fixed at their current best 
values) and the full cycle through all parameters is repeated for a configurable number of outer iterations, with a random perturbation of parameter 
values in between the first and second outer iterations onward.

**************************************************************************************************************************************************************/

// Include the class header
#include "BisectionAlgorithm2.h"


/**************************************************************************************************************************************************************
CTOR

Initializes the algorithm. First, all parameters are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
BisectionAlgorithm::BisectionAlgorithm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read bisection information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Bisection block within the file
        if (CheckToken(pFile, "BeginBisectionAlg", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndBisectionAlg", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginBisectionAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndBisectionAlg") == NULL) {

                if (strstr(line, "ConvergenceVal") != NULL) {
                    // Set the convergence value
                    sscanf(line, "%s %lf", tmp, &m_ConvVal);
                }
                else if (strstr(line, "MaxOuterIterations") != NULL) {
                    // Set the number of outer iterations
                    sscanf(line, "%s %d", tmp, &m_MaxOuter);
                }
                else if (strstr(line, "MaxInnerIterations") != NULL) {
                    // Set the number of inner iterations
                    sscanf(line, "%s %d", tmp, &m_MaxInner);
                }

                // Get the next line
                line = GetNxtDataLine(pFile, inputFilename);
            }
        }
        else {
            // The bisection algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Perform setup checks of the configuration file
    if (m_MaxOuter <= 0) {
        LogError(ERR_FILE_IO, "Invalid number of outer iterations");
        ExitProgram(1);
    }

    if (m_MaxInner <= 0) {
        LogError(ERR_FILE_IO, "Invalid number of inner iterations");
        ExitProgram(1);
    }

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the Bisection algorithm and it's member variables
**************************************************************************************************************************************************************/
void BisectionAlgorithm::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
BuildCandidate()

Builds a full parameter vector by copying the base values and substituting a single parameter's value. Used to construct the independent trial 
evaluations needed at each step of the bisection search.
**************************************************************************************************************************************************************/
std::vector<double> BisectionAlgorithm::BuildCandidate(std::vector<double> baseValues, int paramIndex, double paramValue) {

    std::vector<double> candidate = baseValues;
    candidate[paramIndex] = paramValue;
    return candidate;
}

/**************************************************************************************************************************************************************
WarmStart()

Read the best solution from a previous run.
**************************************************************************************************************************************************************/
void BisectionAlgorithm::WarmStart(void) {
    // todo: create this. Placeholder for restoring the search state (current parameter values, active parameter index, and bracket bounds)
    // from a previously interrupted run.
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using the Bisection algorithm.
**************************************************************************************************************************************************************/
void BisectionAlgorithm::Optimize(void) {

    m_NumParams = m_pParamGroup->GetNumParams();

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    std::vector<int> solveCounts;
    for (int entryWorker = 0; entryWorker < numberOfMpiProcesses; entryWorker++) {
        solveCounts.push_back(0);
    }

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    if (m_bWarmStart) {
        // Load a previous analysis that was interrupted
        std::cout << "Warm start has not been configured for the Bisection algorithm. Exiting the analysis..." << std::endl;
        throw std::invalid_argument("Warm start has not been configured for the Bisection algorithm");
    }

    // Start a clean analysis
    // Write the information to the primary log file
    WriteSetup2("BisectionAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), m_pParamGroup->GetNumParams(), m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();

    // Get the current parameter values, and the bounds of each parameter
    std::vector<double> currentValues;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
        currentValues.push_back(temp->GetInitialValueTransformed());
        lowerBounds.push_back(temp->GetLowerBoundTransformed());
        upperBounds.push_back(temp->GetUpperBoundTransformed());
    }

    // Evaluate the initial point
    std::vector<double> initialObjective = std::vector<double>(1, INFINITY);
    ManageSingleObjectiveIterations(std::vector<std::vector<double>>(1, currentValues), 0, m_NumParams, initialObjective);
    m_AlgCount++;

    double Fcur = initialObjective[0];
    double Fold = Fcur;

    m_BestObjective = Fcur;
    m_BestAlternative = currentValues;

    // Set the values into the pointer group and write the initial configuration to the log file
    for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
        m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(currentValues[entryParameter]);
    }
    WriteRecord2(0, Fcur, Fcur, m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Enter the outer solution loop. One outer iteration corresponds to the application of bisection to every parameter.
    for (int entryOuter = 0; entryOuter < m_MaxOuter; entryOuter++) {

        if (entryOuter > 0) {
            // Randomly perturb every parameter within its bounds ahead of the next cycle through the parameters
            for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
                double range = upperBounds[entryParameter] - lowerBounds[entryParameter];
                double r = (double)MyRand() / (double)MY_RAND_MAX;
                currentValues[entryParameter] = (r * range) + lowerBounds[entryParameter];
            }
        }

        m_CurIter = entryOuter + 1;

        WriteInnerEval2(WRITE_BIS, 5 * m_NumParams + m_MaxInner * m_NumParams * 2, '.');
        int tmp = 0;

        double Fmin = Fcur; // Tracks the minimum objective found while processing the current parameter

        // Process one parameter at a time, holding all other parameters fixed at their current values
        for (int entryParam = 0; entryParam < m_NumParams; entryParam++) {

            double Xcur = currentValues[entryParam];
            double Xupr = upperBounds[entryParam];
            double Xlwr = lowerBounds[entryParam];
            double Xqtr = Xlwr + 0.25 * (Xupr - Xlwr);
            double Xmid = Xlwr + 0.50 * (Xupr - Xlwr);
            double X3qt = Xlwr + 0.75 * (Xupr - Xlwr);
            double Xmin = Xcur;

            // Assign and evaluate the initial points. These subdivide the domain into 4 quadrants. The 5 points are all
            // independent of one another (each is derived purely from the parameter's bounds) so they can be evaluated together.
            std::vector<std::vector<double>> initialCandidates;
            initialCandidates.push_back(BuildCandidate(currentValues, entryParam, Xupr));
            initialCandidates.push_back(BuildCandidate(currentValues, entryParam, Xlwr));
            initialCandidates.push_back(BuildCandidate(currentValues, entryParam, Xqtr));
            initialCandidates.push_back(BuildCandidate(currentValues, entryParam, Xmid));
            initialCandidates.push_back(BuildCandidate(currentValues, entryParam, X3qt));

            std::vector<double> initialObjectives = std::vector<double>(initialCandidates.size(), INFINITY);
            ManageSingleObjectiveIterations(initialCandidates, 0, m_NumParams, initialObjectives);
            m_AlgCount += (int)initialCandidates.size();
            for (int entryEval = 0; entryEval < (int)initialCandidates.size(); entryEval++) { WriteInnerEval2(++tmp, 0, '.'); }

            double Fupr = initialObjectives[0];
            double Flwr = initialObjectives[1];
            double Fqtr = initialObjectives[2];
            double Fmid = initialObjectives[3];
            double F3qt = initialObjectives[4];

            // Perform bisections. Each bisection will reduce the search space by 50%.
            for (int entryInner = 0; entryInner < m_MaxInner; entryInner++) {

                // Track which point(s) need to be freshly evaluated this iteration, along with a label indicating where the
                // result of each evaluation should be assigned once the batch completes.
                std::vector<std::vector<double>> pendingCandidates;
                std::vector<int> pendingLabels; // 0 = mid, 1 = qtr, 2 = 3qt

                // is mid-point the current minimum?
                if ((Fmid <= Fupr) && (Fmid <= Flwr) && (Fmid <= Fqtr) && (Fmid <= F3qt)) {
                    // store new minimum
                    Xmin = Xmid;
                    Fmin = Fmid;

                    // shrink domain
                    Xlwr = Xqtr;
                    Flwr = Fqtr;
                    Xupr = X3qt;
                    Fupr = F3qt;
                }
                // is quarter-point the current minimum?
                else if ((Fqtr <= Fupr) && (Fqtr <= Flwr) && (Fqtr <= Fmid) && (Fqtr <= F3qt)) {
                    // store new minimum
                    Xmin = Xqtr;
                    Fmin = Fqtr;

                    // shrink domain
                    Xupr = Xmid;
                    Fupr = Fmid;
                    Xmid = Xqtr;
                    Fmid = Fqtr;
                }
                // is three-quarter-point the current minimum?
                else if ((F3qt <= Fupr) && (F3qt <= Flwr) && (F3qt <= Fmid) && (F3qt <= Fqtr)) {
                    // store new minimum
                    Xmin = X3qt;
                    Fmin = F3qt;

                    // shrink domain
                    Xlwr = Xmid;
                    Flwr = Fmid;
                    Xmid = X3qt;
                    Fmid = F3qt;
                }
                // is upper-bound the current minimum?
                else if ((Fupr <= F3qt) && (Fupr <= Flwr) && (Fupr <= Fmid) && (Fupr <= Fqtr)) {
                    // store the new minimum
                    Xmin = Xupr;
                    Fmin = Fupr;

                    // shrink domain
                    Xlwr = X3qt;
                    Flwr = F3qt;

                    // the new mid-point does not coincide with any previously evaluated point, so it needs a fresh evaluation
                    Xmid = Xlwr + 0.5 * (Xupr - Xlwr);
                    pendingCandidates.push_back(BuildCandidate(currentValues, entryParam, Xmid));
                    pendingLabels.push_back(0);
                }
                // is lower-bound the current minimum?
                else if ((Flwr <= F3qt) && (Flwr <= Fupr) && (Flwr <= Fmid) && (Flwr <= Fqtr)) {
                    // store the new minimum
                    Xmin = Xlwr;
                    Fmin = Flwr;

                    // shrink domain
                    Xupr = Xqtr;
                    Fupr = Fqtr;

                    // the new mid-point does not coincide with any previously evaluated point, so it needs a fresh evaluation
                    Xmid = Xlwr + 0.5 * (Xupr - Xlwr);
                    pendingCandidates.push_back(BuildCandidate(currentValues, entryParam, Xmid));
                    pendingLabels.push_back(0);
                }
                // assume mid-point
                else {
                    // store new minimum
                    Xmin = Xmid;
                    Fmin = Fmid;

                    // shrink domain
                    Xlwr = Xqtr;
                    Flwr = Fqtr;
                    Xupr = X3qt;
                    Fupr = F3qt;
                }

                // The new quarter points always need to be freshly evaluated, based on the (possibly updated) bounds
                Xqtr = Xlwr + 0.25 * (Xupr - Xlwr);
                pendingCandidates.push_back(BuildCandidate(currentValues, entryParam, Xqtr));
                pendingLabels.push_back(1);

                X3qt = Xlwr + 0.75 * (Xupr - Xlwr);
                pendingCandidates.push_back(BuildCandidate(currentValues, entryParam, X3qt));
                pendingLabels.push_back(2);

                // Evaluate the pending points together, since they are all independent of one another
                std::vector<double> pendingObjectives = std::vector<double>(pendingCandidates.size(), INFINITY);
                ManageSingleObjectiveIterations(pendingCandidates, 0, m_NumParams, pendingObjectives);
                m_AlgCount += (int)pendingCandidates.size();
                for (int entryEval = 0; entryEval < (int)pendingCandidates.size(); entryEval++) { WriteInnerEval2(++tmp, 0, '.'); }

                // Assign the results back to the appropriate point
                for (int entryPending = 0; entryPending < (int)pendingCandidates.size(); entryPending++) {
                    if (pendingLabels[entryPending] == 0) { Fmid = pendingObjectives[entryPending]; }
                    else if (pendingLabels[entryPending] == 1) { Fqtr = pendingObjectives[entryPending]; }
                    else { F3qt = pendingObjectives[entryPending]; }
                }
            }

            if (Fmin < Fcur) {
                Fcur = Fmin;
                currentValues[entryParam] = Xmin;

                // Snapshot the full vector at the moment of acceptance. This keeps the logged/displayed solution
                // (m_BestAlternative) always consistent with the logged objective (m_BestObjective), even though
                // currentValues itself may go on to hold un-accepted, randomly-restarted values for parameters
                // that have not yet been (re-)processed later in this same outer iteration, or in a later outer
                // iteration whose random restart fails to improve on the best solution found so far.
                m_BestObjective = Fcur;
                m_BestAlternative = currentValues;
            }
            // else the search failed for this parameter; currentValues[entryParam] is left unchanged (i.e. reverted to Xcur)
        }
        WriteInnerEval2(WRITE_ENDED, 0, '.');

        // Set the best-known values (not the possibly-still-randomized currentValues) into the pointer group and
        // write the iteration result to the log file
        for (int entryParameter = 0; entryParameter < m_NumParams; entryParameter++) {
            m_pParamGroup->GetParamPtr(entryParameter)->SetEstimatedValueTransformed(m_BestAlternative[entryParameter]);
        }
        WriteRecord2(entryOuter + 1, Fcur, fabs(Fold - Fcur), m_pObsGroup, m_pObjFunc, m_pParamGroup);

        Fold = Fcur;
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
void BisectionAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
    fprintf(pFile, "Max Inner Iterations    : %d\n", m_MaxInner);

    // Write the parameter to the file
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
void BisectionAlgorithm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Bisection Method\n");
    fprintf(pFile, "Max Outer Iters          : %d\n", m_MaxOuter);
    fprintf(pFile, "Max Inner Iters          : %d\n", m_MaxInner);
    fprintf(pFile, "Convergence Val          : %E\n", m_ConvVal);
    fprintf(pFile, "Iterations               : %d\n", m_CurIter);
    fprintf(pFile, "Algorithm Evals          : %d\n", m_AlgCount);

    // Close the log file
    fclose(pFile);
}


/**************************************************************************************************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using the Bisection algorithm.
**************************************************************************************************************************************************************/
/*void BisectionAlgorithm::Calibrate(void) {
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);

   Optimize();

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   if(id == 0) {
      sprintf(fileName, "OstOutput%d.txt", id);

      //write statistics of best parameter set to output file
      pFile = fopen(fileName, "a");
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      //write statistics of best parameter set to output file
      m_pStats->WriteStats(stdout);
   }
}*/

/**************************************************************************************************************************************************************
BIS_Program()

Calibrate the model using the Bisection algorithm.
**************************************************************************************************************************************************************/
void BIS_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the bisection algorithm
       BisectionAlgorithm *BisAlg = new BisectionAlgorithm();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       BisAlg->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
}
