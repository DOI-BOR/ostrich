/**************************************************************************************************************************************************************

Parallel PADDS (Pareto Archived Dynamically Dimensioned Search). This is the parallel version of the PADDS multi-objective algorithm, and shares the PADDS
archive, selection metric and hypervolume functions.

The pre-refactor implementation was asynchronous: each worker was sent a new candidate as soon as it returned a result. The Algorithm base class solves
sets of alternatives using a blocking batch call, so this version is synchronous. Candidates are generated in batches (one per available worker) from the
current archive, solved concurrently, and then added to the archive in the order they were generated. Each candidate is compared against the solution it
was generated from.

**************************************************************************************************************************************************************/

// Include the class header
#include "ParaPADDS2.h"

/**************************************************************************************************************************************************************
CTOR

Reads the parallel PADDS configuration. The first block found is used, in order of precedence: BeginParallelPADDSAlg, BeginParallelPADDS,
BeginParaPADDSAlg, BeginParaPADDS, BeginPADDSAlg, BeginPADDS.
**************************************************************************************************************************************************************/
ParaPADDS::ParaPADDS() : PADDS(false) {

    m_SetupName = "ParaPADDS";
    m_AlgorithmName = "ParaPADDS - Parallel Pareto Archived Dynamically Dimensioned Search";

    ReadConfiguration(std::vector<std::string>{ "BeginParallelPADDSAlg", "BeginParallelPADDS", "BeginParaPADDSAlg", "BeginParaPADDS",
                                                "BeginPADDSAlg", "BeginPADDS" }, false);
}

/**************************************************************************************************************************************************************
Optimize()

Search for the Pareto front using parallel PADDS.
**************************************************************************************************************************************************************/
void ParaPADDS::Optimize(void) {

    // Number of candidates solved together: one per secondary worker, plus one if solving on the primary
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);
    m_BatchSize = (numberOfMpiProcesses - 1) + (m_bSolveOnPrimary ? 1 : 0);
    if (m_BatchSize < 1) m_BatchSize = 1;

    // Validate the setup, configure the workers and write the setup
    BeginAnalysis();

    // Number of initial solutions: at least one per worker, limited to the budget
    int its = (int)(0.005 * m_maxiter);
    if (its < m_BatchSize) its = m_BatchSize;
    its = std::min(its, m_maxiter);

    // Create the random initial solutions, solve them together, then add them to the archive in order
    std::vector<std::vector<double>> samples;
    for (int i = 0; i < its; i++) {
        samples.push_back(RandomSolution());
    }
    std::vector<std::vector<double>> objectives = Evaluate(samples);

    for (size_t i = 0; i < samples.size(); i++) {
        m_CurIter++;
        if (UpdateArchive(samples[i], objectives[i]) == ARCHIVE_NON_DOM) {
            WriteArchiveRecord(m_CurIter, (double)(m_maxiter - m_CurIter));
        }
    }

    // Number of evaluations remaining for the DDS search
    int iLeft = m_maxiter - its;

    // Calculating selection metric Z
    Calc_Z();

    // Main loop
    FILE* pPnFile = fopen("OstPADDSPn.txt", "w");
    fprintf(pPnFile, "EVAL  Pn\n");
    fclose(pPnFile);

    int eval = 0;
    while (eval < iLeft) {
        if (IsQuit() == true) { break; }

        // DDS perturbation probability, based on the number of search evaluations completed. Every decision variable is perturbed for the first evaluations
        // (up to twice the batch size).
        double Pn;
        if (eval <= 2 * m_BatchSize) {
            Pn = 1.00;
        }
        else if (iLeft - 2 * m_BatchSize > 1) {
            Pn = 1.0 - log((double)(eval - 2 * m_BatchSize)) / log((double)(iLeft - 2 * m_BatchSize));
        }
        else {
            Pn = 0.00;
        }

        // Generate a batch of candidates from the current archive
        int batchSize = std::min(m_BatchSize, iLeft - eval);
        std::vector<Solution> parents, candidates;
        std::vector<std::vector<double>> batch;

        pPnFile = fopen("OstPADDSPn.txt", "a");
        for (int k = 0; k < batchSize; k++) {
            Solution sbest = SelectParent();
            Solution stest = PerturbSolution(sbest, Pn);

            parents.push_back(sbest);
            candidates.push_back(stest);
            batch.push_back(stest.X);

            fprintf(pPnFile, "%04d  %E\n", eval + k + 1, Pn);
        }
        fclose(pPnFile);

        // Solve the batch, then add the candidates to the archive in the order they were generated
        objectives = Evaluate(batch);
        for (int k = 0; k < batchSize; k++) {
            candidates[k].F = objectives[k];
            ProcessCandidate(candidates[k], parents[k]);
        }

        eval += batchSize;
    }

    // Terminate the secondary workers and let them clean up.
    TerminateWorkers();

    // Write algorithm metrics and the objectives of the non-dominated solutions
    WriteEndingMetrics();
    WriteFrontObjectives();

} /* end Optimize() */

/**************************************************************************************************************************************************************
WriteFrontObjectives()

Write the objectives of the non-dominated solutions to c_results.txt (tab-separated, one solution per line)
**************************************************************************************************************************************************************/
void ParaPADDS::WriteFrontObjectives(void) {

    FILE* pFile = fopen("c_results.txt", "w");
    for (size_t entry = 0; entry < m_NonDom.size(); entry++) {
        for (int j = 0; j < m_num_objs - 1; j++) {
            fprintf(pFile, "%E\t", m_NonDom[entry].F[j]);
        }
        fprintf(pFile, "%E\n", m_NonDom[entry].F[m_num_objs - 1]);
    }
    fclose(pFile);
}

/**************************************************************************************************************************************************************
PARA_PADDS_Program()

Search for the Pareto front using parallel PADDS.
**************************************************************************************************************************************************************/
void PARA_PADDS_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        ParaPADDS* pParaPADDS = new ParaPADDS();

        // Begin the primary worker operations. This transfers information to the secondary workers and manages the solution process
        pParaPADDS->Optimize();
    }
    else {
        // This is a secondary worker. Create the secondary worker, telling it to configure from MPI, and start the work
        ModelWorker worker = ModelWorker(true);
        worker.WorkMPI();
    }
} /* end PARA_PADDS_Program() */
