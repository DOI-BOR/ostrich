/**************************************************************************************************************************************************************

Particle Swarm Optimization (PSO) applies concepts of social behavior to solve optimization problems. The PSO algorithm starts with a 'swarm' of particles
(solutions) and "flies" this population through the design space in search of the optimal solution. At each iteration, a given particle uses it's own prior
best solution (cognitive behavior) along with the current best solution of all particles (social behavior) to decide where to go next.

**************************************************************************************************************************************************************/

// Include the class header
#include "ParticleSwarm.h"


/**************************************************************************************************************************************************************
CTOR

Initializes the swarm. First, all parameters are assigned default values and then the user input file is checked for overriding values.
**************************************************************************************************************************************************************/
ParticleSwarm::ParticleSwarm() {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();

    // Read particle swarm information from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile != NULL) {
        // Identify the Particle Swarm block within the file
        if (CheckToken(pFile, "BeginParticleSwarm", inputFilename) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndParticleSwarm", inputFilename);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginParticleSwarm", inputFilename);
            line = GetNxtDataLine(pFile, inputFilename);

            while (strstr(line, "EndParticleSwarm") == NULL) {

                if (strstr(line, "SwarmSize") != NULL) {
                    // Set the swarm size
                    sscanf(line, "%s %d", tmp, &m_SwarmSize);
                }
                else if (strstr(line, "NumGenerations") != NULL) {
                    // Set the number of generations
                    sscanf(line, "%s %d", tmp, &m_NumGenerationsMaximum);
                }
                else if (strstr(line, "ConstrictionFactor") != NULL) {
                    // Set the constriction factor
                    sscanf(line, "%s %lf", tmp, &m_Constriction);
                }
                else if (strstr(line, "CognitiveParam") != NULL) {
                    // Set the cognitive weight
                    sscanf(line, "%s %lf", tmp, &m_Cognitive);
                }
                else if (strstr(line, "SocialParam") != NULL) {
                    // Set the social weight
                    sscanf(line, "%s %lf", tmp, &m_Social);
                }
                else if (strstr(line, "InertiaWeight") != NULL) {
                    // Set the inertia weight
                    sscanf(line, "%s %lf", tmp, &m_Inertia);
                }
                else if (strstr(line, "InertiaReductionRate") != NULL) {
                    // Set the inertia reduction rate, or flag linear reduction to zero
                    sscanf(line, "%s %s", tmp, tmp2);
                    MyStrLwr(tmp2);
                    if (strcmp(tmp2, "linear") == 0) { m_bLinearInertiaReduction = true; }
                    else { sscanf(line, "%s %lf", tmp, &m_InertiaReductionRate); }
                }
                else if (strstr(line, "InitPopulationMethod") != NULL) {
                    // Set the initializaiton method
                    sscanf(line, "%s %s", tmp, tmp2);
                    MyStrLwr(tmp2);
                    if (strcmp(tmp2, "random") == 0) { m_InitType = RANDOM_INIT; }
                    else if (strcmp(tmp2, "quadtree") == 0) { m_InitType = QUAD_TREE_INIT; }
                    else if (strcmp(tmp2, "lhs") == 0) { m_InitType = LHS_INIT; }
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
            // The particle swarm algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        fclose(pFile);
    }

    // Perform setup checks of the configuration file
    // Check that the swarm size is greater than zero
    if (m_SwarmSize <= 0) {
        LogError(ERR_FILE_IO, "Invalid swarm size");
        ExitProgram(1);
    }

    // Check that the number of generations is greater than zero
    if (m_NumGenerationsMaximum <= 0) {
        LogError(ERR_FILE_IO, "Invalid number of generations");
        ExitProgram(1);
    }

    // Remember the starting inertia weight so that it can be reduced each generation and restored afterwards
    m_InertiaInitial = m_Inertia;

    // TODO: port support for seeding specific swarm members via a BeginInitParams/EndInitParams block,
    // similar to what the pre-refactor ParticleSwarm supported. This has not yet been ported to the
    // new Algorithm-based structure (GeneticAlgorithm has the same limitation currently).

    IncCtorCount();
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by the PSO and it's member variables
**************************************************************************************************************************************************************/
void ParticleSwarm::Destroy(void) {

    delete [] m_pTrees;

    m_Positions.clear();
    m_Velocities.clear();
    m_PersonalBestPositions.clear();
    m_PersonalBestObjectives.clear();

    IncDtorCount();
}

/**************************************************************************************************************************************************************
CalcSwarmMedian()

Calculates and returns the median objective function of the swarm. This value is used in the termination criteria of the PSO algorithm.
**************************************************************************************************************************************************************/
double ParticleSwarm::CalcSwarmMedian(std::vector<double> objectives) {

    return CalcMedian(objectives.data(), (int)objectives.size());
}

/**************************************************************************************************************************************************************
CreateInitialSwarm()

Creates an initial set of particle positions for use in the particle swarm algorithm.
**************************************************************************************************************************************************************/
std::vector<std::vector<double>> ParticleSwarm::CreateInitialSwarm(int swarmSize) {

    // Create the holder for the samples
    int numberOfParameters = m_pParamGroup->GetNumParams();
    std::vector<std::vector<double>> samples;

    // Generate the set of initial particle positions
    if (m_InitType == LHS_INIT) {
        // Use Latin Hypercube to initialize the space
        LatinHypercube* pLHS = new LatinHypercube(swarmSize, numberOfParameters, numberOfParameters);
        MEM_CHECK(pLHS);

        // Extract the bounds needed to create the LHS matrix
        double* lowerBounds = new double[numberOfParameters];
        double* upperBounds = new double[numberOfParameters];

        for (int i = 0; i < numberOfParameters; i++) {
            lowerBounds[i] = m_pParamGroup->GetParamPtr(i)->GetLowerBoundTransformed();
            upperBounds[i] = m_pParamGroup->GetParamPtr(i)->GetUpperBoundTransformed();
        }

        // Fill the array from the LHS object
        pLHS->CreateUniformSample(lowerBounds, upperBounds);
        samples = pLHS->GetSampleMatrix();

        // Cleanup the memory
        delete pLHS;
        delete[] lowerBounds;
        delete[] upperBounds;

    } else if (m_InitType == QUAD_TREE_INIT) {
        // Use quad trees to initialize the space
        // Initialize the trees if needed
        if (m_pTrees == NULL) {
            m_TreeSize = numberOfParameters;
            m_pTrees = new QuadTree[m_TreeSize];
            MEM_CHECK(m_pTrees);
            for (int j = 0; j < m_TreeSize; j++) {
                double lwr = m_pParamGroup->GetParamPtr(j)->GetLowerBoundTransformed();
                double upr = m_pParamGroup->GetParamPtr(j)->GetUpperBoundTransformed();
                m_pTrees[j].Init(lwr, upr);
            }
        }

        int lvl = 0;
        int idx = 0;
        for (int entryParticle = 0; entryParticle < swarmSize; entryParticle++) {
            double* pVals = GetTreeCombo(lvl, idx, m_pTrees, m_TreeSize);

            // Expand the tree if needed
            if (pVals == NULL) {
                for (int j = 0; j < m_TreeSize; j++) { m_pTrees[j].Expand(); }
                lvl++;
                idx = 0;
                pVals = GetTreeCombo(lvl, idx, m_pTrees, m_TreeSize);
            }
            idx++;

            // Copy the tree values into a chromosome for the swarm
            std::vector<double> temp;
            for (int j = 0; j < numberOfParameters; j++) {
                temp.push_back(pVals[j]);
            }
            samples.push_back(temp);

            delete[] pVals;
        }

    } else if (m_InitType == RANDOM_INIT) {
        // Use random values to initialize the space
        for (int entryParticle = 0; entryParticle < swarmSize; entryParticle++) {
            // Create the vector for this particle
            std::vector<double> temp;

            // Loop and fill the vector
            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                double lwr = m_pParamGroup->GetParamPtr(entryParameter)->GetLowerBoundTransformed();
                double upr = m_pParamGroup->GetParamPtr(entryParameter)->GetUpperBoundTransformed();
                double range = upr - lwr;
                double r = (double)MyRand() / (double)MY_RAND_MAX;

                temp.push_back((r * range) + lwr);
            }

            // Append the vector to the sample vector
            samples.push_back(temp);
        }
    }
    else {
        // Log an error
        LogError(ERR_FILE_IO, "Unknown particle swarm initialization method");
        ExitProgram(1);
    }

    // Return to the calling function
    return samples;
}

/**************************************************************************************************************************************************************
ConstrainToBounds()

Constrains a particle's revised position to stay within its parameter bounds, taking care to preserve the angle (i.e. direction of movement) of the particle
by scaling back the velocity of every parameter of the particle by the same fraction.
**************************************************************************************************************************************************************/
void ParticleSwarm::ConstrainToBounds(int particleIndex) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Determine the most restrictive fractional move needed to keep every parameter of the particle in bounds
    double dxMin = 1.00;
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        double lwr = m_pParamGroup->GetParamPtr(entryParameter)->GetLowerBoundTransformed();
        double upr = m_pParamGroup->GetParamPtr(entryParameter)->GetUpperBoundTransformed();

        double v = m_Velocities[particleIndex][entryParameter];
        double x = m_Positions[particleIndex][entryParameter] - v; // original (pre-move) position

        if (m_Positions[particleIndex][entryParameter] > upr) {
            double dxOld = v;
            double dxNew = 0.5 * (upr - x);
            double dxFrac = fabs(dxNew / dxOld); // relative change
            if (dxFrac < dxMin) { dxMin = dxFrac; }
            m_NumUprViols++;
        }
        if (m_Positions[particleIndex][entryParameter] < lwr) {
            double dxOld = v;
            double dxNew = 0.5 * (lwr - x);
            double dxFrac = fabs(dxNew / dxOld); // relative change
            if (dxFrac < dxMin) { dxMin = dxFrac; }
            m_NumLwrViols++;
        }
    }

    // Apply the most restrictive fraction to every parameter of the particle
    for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
        double v = m_Velocities[particleIndex][entryParameter];
        double x = m_Positions[particleIndex][entryParameter] - v; // original (pre-move) position

        m_Velocities[particleIndex][entryParameter] *= dxMin; // revised velocity
        m_Positions[particleIndex][entryParameter] = x + (v * dxMin); // revised position
    }
}

/**************************************************************************************************************************************************************
UpdateSwarm()

Updates the velocity and position of each particle based on it's own personal best solution (cognitive behavior) and the swarm's current best solution
(social behavior), then constrains each particle to remain within the parameter bounds.
**************************************************************************************************************************************************************/
void ParticleSwarm::UpdateSwarm(int generation) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    for (int entryParticle = 0; entryParticle < m_SwarmSize; entryParticle++) {
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            // Intermediary variables
            double x = m_Positions[entryParticle][entryParameter];
            double pl = m_PersonalBestPositions[entryParticle][entryParameter];
            double pg = m_BestAlternative[entryParameter];

            // Random weights
            double r1 = (double)MyRand() / (double)MY_RAND_MAX;
            double r2 = (double)MyRand() / (double)MY_RAND_MAX;

            // Revised velocity
            double v = m_Velocities[entryParticle][entryParameter];
            v = m_Constriction * ((m_Inertia * v) + m_Cognitive * r1 * (pl - x) + m_Social * r2 * (pg - x));

            // Assign minimum perturbation to prevent stagnation
            double vmin;
            if (strcmp(m_pParamGroup->GetParamPtr(entryParameter)->GetType(), "real") == 0) {
                vmin = (0.01 * fabs(x)) / (generation + 1);
            } else {
                vmin = 0.50;
            }

            if (fabs(v) < vmin) {
                // Adjust randomized minimum velocity
                double sgn = (double)MyRand() / (double)MY_RAND_MAX; // random direction
                if (sgn >= 0.50) { v = +((1.00 + r1) * vmin); }
                else { v = -((1.00 + r2) * vmin); }
            }
            m_Velocities[entryParticle][entryParameter] = v;

            // Revised position
            m_Positions[entryParticle][entryParameter] = x + v;
        }

        // Constrain the revised position to stay within parameter limits, preserving direction of movement
        ConstrainToBounds(entryParticle);
    }

    // Reduce the inertia weight for the next generation
    if (m_bLinearInertiaReduction) {
        // Linearly reduce to zero over the course of the run
        m_Inertia = m_InertiaInitial;
        m_InertiaReductionRate = (double)generation / (double)m_NumGenerationsMaximum;
    }
    m_Inertia *= (1.00 - m_InertiaReductionRate);
}

/**************************************************************************************************************************************************************
WarmStart()

Read the best solution from a previous run.
**************************************************************************************************************************************************************/
void ParticleSwarm::WarmStart(void) {
    // todo: create this. Placeholder for restoring swarm positions/velocities and the global/personal best
    // solutions from a previously interrupted run.
}

/**************************************************************************************************************************************************************
Optimize()

Minimize the objective function using PSO.
**************************************************************************************************************************************************************/
void ParticleSwarm::Optimize(void) {

    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Fill the solve counter with zeros
    std::vector<int> solveCounts = std::vector<int>(numberOfMpiProcesses, 0);

    // Initialize the workers
    ConfigureWorkers(solveCounts);

    if (m_bWarmStart) {
        // Load a previous analysis that was interrupted
        std::cout << "Warm start has not been configured for the Particle Swarm algorithm. Exiting the analysis..." << std::endl;
        throw std::invalid_argument("Warm start has not been configured for the Particle Swarm algorithm");

    } else {
        // Start a clean analysis
        // Write the information to the primary log file
        WriteSetup2("ParticleSwarm", m_ExecCmd, m_pObjFunc->GetObjFuncStr(), m_pParamGroup->GetNumParams(), m_pParamGroup->GetNumTiedParams());
        WriteStartingMetrics();

        // Generate the remaining swarm positions
        m_Positions = CreateInitialSwarm(m_SwarmSize - 1);

        // Initialize the current parameter values as the first particle
        std::vector<double> initialConditions;
        for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParameter);
            initialConditions.push_back(temp->GetInitialValueTransformed());
        }
        m_Positions.insert(m_Positions.begin(), initialConditions);
    }

    // Initial velocities are zero
    m_Velocities = std::vector<std::vector<double>>(m_SwarmSize, std::vector<double>(numberOfParameters, 0.00));

    // Initialize the objectives for the first solve
    std::vector<double> objectives = std::vector<double>(m_SwarmSize, INFINITY);

    // Evaluate the initial swarm
    ManageSingleObjectiveIterations(m_Positions, 0, numberOfParameters, objectives);

    // Personal bests start out equal to the initial positions/objectives
    m_PersonalBestPositions = m_Positions;
    m_PersonalBestObjectives = objectives;

    // Compute the swarm's best and set it as the running best, if it is an improvement
    GetBestSingleObjective(objectives, m_BestObjectiveIteration, m_BestObjectiveIndexIteration);
    if (m_BestObjectiveIteration < m_BestObjective) {
        m_BestObjective = m_BestObjectiveIteration;
        m_BestAlternative = m_Positions[m_BestObjectiveIndexIteration];
    }

    double median = CalcSwarmMedian(objectives);
    m_CurStop = fabs((median - m_BestObjective) / median);

    // Set the values into the pointer groups
    for (int entryParam = 0; entryParam < numberOfParameters; entryParam++) {
        ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParam);
        temp->SetEstimatedValueTransformed(m_BestAlternative[entryParam]);
    }

    // Write the initial configuration to the log file
    WriteRecord2(m_Generation, m_BestObjective, m_CurStop, m_pObsGroup, m_pObjFunc, m_pParamGroup);

    // Enter the solution loop
    while (m_NumGenerationsMaximum > m_Generation && m_CurStop >= m_StopVal) {

        // Update the velocity and position of every particle
        UpdateSwarm(m_Generation);

        // Evaluate the updated swarm
        ManageSingleObjectiveIterations(m_Positions, 0, numberOfParameters, objectives);

        // Revise personal and global bests
        for (int entryParticle = 0; entryParticle < m_SwarmSize; entryParticle++) {
            if (objectives[entryParticle] < m_PersonalBestObjectives[entryParticle]) {
                m_PersonalBestPositions[entryParticle] = m_Positions[entryParticle];
                m_PersonalBestObjectives[entryParticle] = objectives[entryParticle];
            }
        }

        GetBestSingleObjective(m_PersonalBestObjectives, m_BestObjectiveIteration, m_BestObjectiveIndexIteration);
        if (m_BestObjectiveIteration < m_BestObjective) {
            m_BestObjective = m_BestObjectiveIteration;
            m_BestAlternative = m_PersonalBestPositions[m_BestObjectiveIndexIteration];
        }

        median = CalcSwarmMedian(objectives);
        m_CurStop = fabs((median - m_BestObjective) / median);

        // Set the values into the pointer groups
        for (int entryParam = 0; entryParam < numberOfParameters; entryParam++) {
            ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParam);
            temp->SetEstimatedValueTransformed(m_BestAlternative[entryParam]);
        }

        // Increment the generation variable
        m_Generation++;

        // Write the iteration result to the log file
        WriteRecord2(m_Generation, m_BestObjective, m_CurStop, m_pObsGroup, m_pObjFunc, m_pParamGroup);
    }

    // Restore the initial inertia weight
    m_Inertia = m_InertiaInitial;

    // Terminate the secondary workers and let them clean up
    TerminateWorkers();

    // Write algorithm metrics
    WriteEndingMetrics();

}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out setup and metrics for the swarm.
**************************************************************************************************************************************************************/
void ParticleSwarm::WriteStartingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Swarm Size              : %d\n", m_SwarmSize);
    fprintf(pFile, "Constriction Factor     : %.2lf\n", m_Constriction);
    fprintf(pFile, "Cognitive Weight        : %.2lf\n", m_Cognitive);
    fprintf(pFile, "Social Weight           : %.2lf\n", m_Social);
    fprintf(pFile, "Inertia Weight          : %.2lf\n", m_Inertia);
    fprintf(pFile, "Initialization Method   : ");
    if (m_InitType == RANDOM_INIT) { fprintf(pFile, "Random\n"); }
    else if (m_InitType == QUAD_TREE_INIT) { fprintf(pFile, "Quad-Tree\n"); }
    else if (m_InitType == LHS_INIT) { fprintf(pFile, "Latin Hypercube Sampling\n"); }
    else { fprintf(pFile, "Unknown\n\n"); }

    // Write the parameter to the file
    fprintf(pFile, "Run   "); // TODO: align columns
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
void ParticleSwarm::WriteEndingMetrics(void) {

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);

    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : Particle Swarm Optimization\n");
    fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
    fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
    fprintf(pFile, "Max Generations         : %d\n", m_NumGenerationsMaximum);
    fprintf(pFile, "Actual Generations      : %d\n", m_Generation);
    fprintf(pFile, "Upper Bound Violations  : %d\n", m_NumUprViols);
    fprintf(pFile, "Lower Bound Violations  : %d\n", m_NumLwrViols);

    if (m_CurStop <= m_StopVal) {
        fprintf(pFile, "Algorithm successfully converged on a solution\n");
    }
    else {
        fprintf(pFile, "Algorithm failed to converge on a solution, more generations may be needed\n");
    }

    // Close the log file
    fclose(pFile);
}


/**************************************************************************************************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PSO.
**************************************************************************************************************************************************************/
/*void ParticleSwarm::Calibrate(void) {
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
PSO_Program()

Calibrate the model using PSO.
**************************************************************************************************************************************************************/
void PSO_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the particle swarm algorithm
       ParticleSwarm *PSO = new ParticleSwarm();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       PSO->Optimize();

   } else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker, telling it to configure from MPI
       ModelWorker worker = ModelWorker(true);

       // Start the work in the processor
       worker.WorkMPI();

       // Tear-down is triggered by the primary worker. Perform any additional cleanup actions
   }
}

