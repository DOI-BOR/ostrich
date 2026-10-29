/**************************************************************************************************************************************************************

This algorithm will "grid" the objective function surface by evaluating an exhaustive (or nearly exhaustive) set of parameter configurations. Parameter
sets are generated deterministically across a user-specified (or default) number of grid points per parameter, and are dispatched for evaluation in
batches ("mini-grids") that are solved using the shared, asynchronous MPI work-farming functionality provided by the Algorithm base class.

**************************************************************************************************************************************************************/

// Include the class header
#include "GridAlgorithm2.h"

/**************************************************************************************************************************************************************
CTOR

Initializes the grid. First, default grid dimensions (two points per parameter) are assigned, and then the user input file is checked for overriding
values. The grid geometry (per-parameter lower bound, step size and rollover counts used to decode a flat grid index) is computed from the parameter
bounds once the grid dimensions are known.
**************************************************************************************************************************************************************/
GridAlgorithm::GridAlgorithm() {

    FILE* pFile;
    char* pTok;
    char* line;
    char tmp[DEF_STR_SZ];
    int numberOfParameters = m_pParamGroup->GetNumParams();
    IroncladString inputFilename = GetInFileName();

    // Assign the default grid dimensions (two points per parameter)
    m_GridDims = std::vector<int>(numberOfParameters, 2);

    // Read grid algorithm information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Grid Algorithm block within the file
        if (CheckToken(pFile, "BeginGridAlg", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndGridAlg", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginGridAlg", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndGridAlg") == NULL) {

                if (strstr(line, "Dimensions") != NULL) {
                    // Extract the number of grid points for each parameter, one-by-one
                    pTok = &line[10];
                    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                        int len = ExtractString(pTok, tmp);
                        len = ValidateExtraction(len, entryParameter, numberOfParameters, "GridAlgorithm::CTOR()");
                        pTok += len;

                        m_GridDims[entryParameter] = atoi(tmp);

                        if (m_GridDims[entryParameter] < 2) {
                            LogError(ERR_FILE_IO, "Invalid grid dimension, defaulting to 2");
                            m_GridDims[entryParameter] = 2;
                        }
                    }
                }
                else if (strstr(line, "EvalsPerIter") != NULL) {
                    // Set the number of grid points evaluated per batch/iteration
                    long long evalsPerIter = 0;
                    sscanf(line, "%s %lld", tmp, &evalsPerIter);
                    m_MiniSize = evalsPerIter;
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
            // The grid algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Determine the total grid size and check validity
    m_GridSize = 1;
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        m_GridSize *= m_GridDims[entryParameter];
    }

    if (m_GridSize < 1) {
        LogError(ERR_FILE_IO, "Invalid grid size");
        ExitProgram(1);
    }

    // Assign a default batch size if one was not specified (or was invalid)
    if (m_MiniSize <= 0) {
        int numberOfMpiProcesses;
        MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);
        m_MiniSize = (long long)10 * numberOfParameters * numberOfMpiProcesses;
    }
    
    if (m_MiniSize > m_GridSize) { m_MiniSize = m_GridSize; }

    // Determine the number of batches/iterations required to cover the entire grid
    m_NumIters = (int)(m_GridSize / m_MiniSize);
    if ((m_GridSize % m_MiniSize) != 0) { m_NumIters++; }

    // Compute the per-parameter lower bound and step size, and the rollover counts used to
    // decode a flat grid index into per-parameter indices (mixed-radix decomposition)
    m_Lwr = std::vector<double>(numberOfParameters, 0.00);
    m_StepSize = std::vector<double>(numberOfParameters, 0.00);
    m_Rval = std::vector<long long>(numberOfParameters, 0);

    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
        double upr = temp->GetUpperBoundTransformed();
        double lwr = temp->GetLowerBoundTransformed();
        double dim = (double)m_GridDims[entryParameter];

        m_Lwr[entryParameter] = lwr;
        m_StepSize[entryParameter] = (upr - lwr) / (dim - 1.00);
    }

    m_Rval[numberOfParameters - 1] = m_GridDims[numberOfParameters - 1];
    for (int entryParameter = numberOfParameters - 2; entryParameter >= 0; entryParameter--) {
        m_Rval[entryParameter] = m_Rval[entryParameter + 1] * m_GridDims[entryParameter];
    }

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the Grid algorithm and it's member variables
**************************************************************************************************************************************************************/
void GridAlgorithm::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
WarmStart()

Read the best solution from a previous run. Not yet supported for the Grid algorithm.
**************************************************************************************************************************************************************/
void GridAlgorithm::WarmStart(void) {
    std::cout << "Warm start has not been configured for the Grid Algorithm. Exiting the analysis..." << std::endl;
    throw std::invalid_argument("Warm start has not been configured for the Grid Algorithm");
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using the Grid algorithm (exhaustive search). Grid points are generated and evaluated in batches ("mini-grids"), with
each batch dispatched to the available MPI workers using the shared work-farming functionality provided by the Algorithm base class.
**************************************************************************************************************************************************************/
void GridAlgorithm::Optimize(void) {

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Warm start is not yet supported for the Grid algorithm
    if (m_bWarmStart) {
        WarmStart();
    }

    // Fill the solve counter with zeros (no warm start support, yet)
    std::vector<int> solveCounts;
    for (int entryWorker = 0; entryWorker < numberOfMpiProcesses; entryWorker++) {
        solveCounts.push_back(0);
    }

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Write the information to the primary log file
    WriteSetup2("GridAlgorithm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), numberOfParameters, m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();

    // Main optimization loop. Each iteration evaluates a batch ("mini-grid") of grid points.
    long long numLeft = m_GridSize;
    long long offset = 0;
    m_CurIter = 0;

    while (numLeft > 0) {
        m_CurIter++;

        long long batchSize = (numLeft < m_MiniSize) ? numLeft : m_MiniSize;

        // Generate the batch of grid points to be solved
        std::vector<std::vector<double>> samples;
        for (long long entryPoint = 0; entryPoint < batchSize; entryPoint++) {
            long long gridIndex = offset + entryPoint;
            std::vector<double> point;

            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                point.push_back(GetGridVal(gridIndex, entryParameter));
            }

            samples.push_back(point);
        }

        // Evaluate the batch. This is dispatched across all available MPI workers by the base class.
        std::vector<double> objectives((size_t)batchSize, INFINITY);
        ManageSingleObjectiveIterations(samples, 0, numberOfParameters, objectives);

        // Determine the best and median objectives within the batch
        double batchBestObjective;
        int batchBestIndex;
        GetBestSingleObjective(objectives, batchBestObjective, batchBestIndex);

        // CalcMedian() sorts its input in place, so operate on a scratch copy
        std::vector<double> objectivesScratch = objectives;
        double median = CalcMedian(objectivesScratch.data(), (int)batchSize);

        // Update the global best, if improved by this batch
        if (batchBestObjective < m_BestObjective) {
            m_BestObjective = batchBestObjective;
            m_BestAlternative = samples[batchBestIndex];
        }

        // Set the current best values into the parameter group so they are reflected in the log/model state
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
            temp->SetEstimatedValueTransformed(m_BestAlternative[entryParameter]);
        }

        // Write the iteration result to the log file
        WriteRecord2(m_CurIter, m_BestObjective, median, m_pObsGroup, m_pObjFunc, m_pParamGroup);

        // Advance to the next batch
        numLeft -= batchSize;
        offset += batchSize;
    }

    // Terminate the secondary workers and let them clean up.
    TerminateWorkers();

    // Write algorithm metrics
    WriteEndingMetrics();

} /* end Optimize() */

/**************************************************************************************************************************************************************
GetGridVal()

Compute value of parameter (paramIndex) at flat grid location (gridIndex), via mixed-radix decomposition of the flat index using the rollover counts
computed from the grid dimensions.
**************************************************************************************************************************************************************/
double GridAlgorithm::GetGridVal(long long gridIndex, int paramIndex) {
    int last = (int)m_GridDims.size() - 1;
    long long idx;

    if (paramIndex == last) {
        idx = gridIndex % m_Rval[paramIndex];
    }
    else {
        idx = (gridIndex % m_Rval[paramIndex]) / m_Rval[paramIndex + 1];
    }

    return (m_Lwr[paramIndex] + m_StepSize[paramIndex] * (double)idx);
} /* end GetGridVal() */

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for the grid.
**************************************************************************************************************************************************************/
void GridAlgorithm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Log the algorithm setup information
    fprintf(pFile, "Grid Size               : %lld\n", m_GridSize);

    fprintf(pFile, "Grid Dimensions         : ");
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        fprintf(pFile, "%d", m_GridDims[entryParameter]);
        if (entryParameter < (numberOfParameters - 1)) fprintf(pFile, " by ");
    }
    fprintf(pFile, "\n");

    fprintf(pFile, "Evals Per Iteration     : %lld\n", m_MiniSize);
    fprintf(pFile, "Max Iterations          : %d\n", m_NumIters);

    // Write the parameter header to the file
    fprintf(pFile, "Run   ");
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Median\n");

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
**************************************************************************************************************************************************************/
void GridAlgorithm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Grid Algorithm (Exhaustive Search)\n");
    fprintf(pFile, "Max Iterations          : %d\n", m_NumIters);
    fprintf(pFile, "Actual Iterations       : %d\n", m_CurIter);
    fprintf(pFile, "Grid Size               : %lld\n", m_GridSize);
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);

    // Close the log file
    fclose(pFile);
}

/**************************************************************************************************************************************************************
GRID_Program()

Optimize the model using the Grid algorithm (exhaustive search).
**************************************************************************************************************************************************************/
void GRID_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the grid algorithm
       GridAlgorithm *Grid = new GridAlgorithm();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       Grid->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
} /* end GRID_Program() */
