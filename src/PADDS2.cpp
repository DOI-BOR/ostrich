/**************************************************************************************************************************************************************

The PADDS (Pareto Archived Dynamically Dimensioned Search) algorithm is a multi-objective version of the DDS algorithm (Asadzadeh and Tolson, 2009, 2013).
It was ported to OSTRICH by L. Shawn Matott from a C++ implementation of Mohammadamin Jahanpour (mjahanpo@uwaterloo.ca).

PADDS maintains an archive of non-dominated and dominated solutions. The initial solutions (user-specified and random) are solved together using the MPI
work-farming functionality provided by the Algorithm base class. Each subsequent candidate is generated from a member of the non-dominated archive and
depends on the result of the previous candidate, so the main search solves one candidate at a time. The objectives are the WSSE of each observation group.

---------------------------------------------------------------------------------------------------------------------------------------------------------------
The copyright, disclaimer, and citation given below affects the following list of subroutines:

   bool covers(double* cub, double * regLow);
   bool partCovers(double* cub, double * regUp);
   int containsBoundary(double* cub, double * regLow, int split);
   double getMeasure(double * regLow, double * regUp);
   int isPile(double* cub, double * regLow, double * regUp);
   double getMedian(double * bounds, int size);
   double computeTrellis(double * regLow, double * regUp, double * trellis);
   void stream(double * regionLow, double * regionUp, double ** points, int npoints, int split, double cover);

Copyright (c) 2006 Nicola Beume  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version. This program is distributed in the
hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.

This program calculates the dominated hypervolume or S-metric of a set of d-dimensional points (d>=3). Please refer to the following publication for a
description of the algorithm:
   Nicola Beume and Guenter Rudolph. Faster S-Metric Calculation by Considering Dominated Hypervolume as Klee's Measure Problem. In: B. Kovalerchuk (ed.):
   Proceedings of the Second IASTED Conference on Computational Intelligence (CI 2006), pp. 231-236. ACTA Press: Anaheim, 2006.

   Extended version published as: Technical Report of the Collaborative Research Centre 531 'Computational Intelligence', CI-216/06, ISSN 1433-3325.
   University of Dortmund, July 2006.
---------------------------------------------------------------------------------------------------------------------------------------------------------------

**************************************************************************************************************************************************************/

// Include the class header
#include "PADDS2.h"
#include "WriteUtility2.h"

/**************************************************************************************************************************************************************
CTOR

Reads the PADDS configuration from the BeginPADDSAlg (or BeginPADDS) block of the OSTRICH input file, and any initial parameter sets from the
BeginInitParams block. Default values are used for any settings that are not specified.
**************************************************************************************************************************************************************/
PADDS::PADDS() {

    ReadConfiguration(std::vector<std::string>{ "BeginPADDSAlg", "BeginPADDS" }, true);

    IncCtorCount();
}

/**************************************************************************************************************************************************************
CTOR

Constructor for derived algorithms (e.g. parallel PADDS), which read their own configuration using ReadConfiguration().
**************************************************************************************************************************************************************/
PADDS::PADDS(bool bReadConfiguration) {

    if (bReadConfiguration) {
        ReadConfiguration(std::vector<std::string>{ "BeginPADDSAlg", "BeginPADDS" }, true);
    }

    IncCtorCount();
}

/**************************************************************************************************************************************************************
ReadConfiguration()

Reads the algorithm configuration from the first of the given blocks that is present in the OSTRICH input file (e.g. "BeginPADDSAlg"; the matching end
tag is formed by replacing "Begin" with "End"). If none are present, the last block name is used and the defaults are kept. Initial parameter sets are
read from the BeginInitParams block if requested.
**************************************************************************************************************************************************************/
void PADDS::ReadConfiguration(std::vector<std::string> blockNames, bool bReadInitParams) {

    FILE* pFile;
    char* line;
    char tmp[DEF_STR_SZ];
    char begin_token[DEF_STR_SZ];
    char end_token[DEF_STR_SZ];
    IroncladString inputFilename = GetInFileName();
    int numberOfParameters = m_pParamGroup->GetNumParams();

    // Read the configuration from the OSTRICH input file
    pFile = fopen(inputFilename, "r");
    if (pFile == NULL) {
        // Couldn't open the file. Use the defaults and log the error.
        sprintf(tmp, "Couldn't open %s config. file. Using Defaults", m_SetupName.c_str());
        LogError(ERR_FILE_IO, tmp);
        return;
    }

    // Accept multiple section headings, in order of precedence
    strcpy(begin_token, blockNames.back().c_str());
    for (size_t entry = 0; entry < blockNames.size(); entry++) {
        rewind(pFile);
        if (CheckToken(pFile, blockNames[entry].c_str(), inputFilename) == true) {
            strcpy(begin_token, blockNames[entry].c_str());
            break;
        }
    }
    sprintf(end_token, "End%s", &begin_token[5]);
    rewind(pFile);

    // Identify the algorithm block within the file
    if (CheckToken(pFile, begin_token, inputFilename) == true) {
        // Identify the end of the algorithm parameter block
        FindToken(pFile, end_token, inputFilename);
        rewind(pFile);

        // Import data starting from the beginning of the block
        FindToken(pFile, begin_token, inputFilename);
        line = GetNxtDataLine(pFile, inputFilename);

        while (strstr(line, end_token) == NULL) {
            if (strstr(line, "PerturbationValue") != NULL) {
                sscanf(line, "%s %lf", tmp, &m_fraction1);
            }
            else if (strstr(line, "MaxIterations") != NULL) {
                sscanf(line, "%s %d", tmp, &m_maxiter);
            }
            else if (strstr(line, "SelectionMetric") != NULL) {
                sscanf(line, "SelectionMetric %s", tmp);

                MyStrLwr(tmp);
                if (strcmp(tmp, "random") == 0) m_Select_metric = 0;
                else if (strcmp(tmp, "crowdingdistance") == 0) m_Select_metric = 1;
                else if (strcmp(tmp, "estimatedhypervolumecontribution") == 0) m_Select_metric = 2;
                else if (strcmp(tmp, "exacthypervolumecontribution") == 0) m_Select_metric = 3;
                else m_Select_metric = 3;
            }

            // Get the next line
            line = GetNxtDataLine(pFile, inputFilename);
        }
    }
    else {
        // The algorithm has been specified but there is no input block. Log an error.
        sprintf(tmp, "Using default %s setup.", m_SetupName.c_str());
        LogError(ERR_FILE_IO, tmp);
    }

    // Read in a list of initial parameter sets
    rewind(pFile);
    if (bReadInitParams && CheckToken(pFile, "BeginInitParams", inputFilename) == true) {
        FindToken(pFile, "EndInitParams", inputFilename);
        rewind(pFile);

        FindToken(pFile, "BeginInitParams", inputFilename);
        line = GetNxtDataLine(pFile, inputFilename);
        while (strstr(line, "EndInitParams") == NULL) {
            // Extract the values, one-by-one, making any necessary conversions
            char* pTok = line;
            std::vector<double> initial;
            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                int j = ExtractString(pTok, tmp);
                j = ValidateExtraction(j, entryParameter, numberOfParameters, "PADDS::ReadConfiguration()");
                pTok += j;
                initial.push_back(m_pParamGroup->GetParamPtr(entryParameter)->ConvertInVal(atof(tmp)));
            }
            m_Init.push_back(initial);

            line = GetNxtDataLine(pFile, inputFilename);
        }
    }

    fclose(pFile);
}

/**************************************************************************************************************************************************************
Destroy()

Free up memory used by PADDS and it's member variables
**************************************************************************************************************************************************************/
void PADDS::Destroy(void) {

    IncDtorCount();
}

/**************************************************************************************************************************************************************
BeginAnalysis()

Validate the configuration, configure the multi-objective solution and the workers, and write the setup to the primary log file
**************************************************************************************************************************************************************/
void PADDS::BeginAnalysis(void) {

    // Get MPI information
    int numberOfMpiProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfMpiProcesses);

    // Warm start is not available for PADDS. Log the request and continue with a new analysis.
    if (m_bWarmStart) {
        char msg[DEF_STR_SZ];
        sprintf(msg, "Warm start is not available for the %s algorithm. A new analysis will be performed.", m_SetupName.c_str());
        LogError(ERR_FILE_IO, msg);
        std::cout << msg << std::endl;
    }

    // The workers compute the objectives as the WSSE of each observation group. Other objective functions are not yet supported.
    if (m_ObjFuncId != OBJ_FUNC_WSSE) {
        char msg[DEF_STR_SZ];
        sprintf(msg, "%s requires the WSSE objective function (one objective per observation group).", m_SetupName.c_str());
        LogError(ERR_FILE_IO, msg);
        std::cout << msg << " Exiting the analysis..." << std::endl;
        ExitProgram(1);
    }

    // Validate the configuration
    if (m_maxiter < 1) {
        LogError(ERR_FILE_IO, "Invalid PADDS MaxIterations. Defaulting to 50.");
        m_maxiter = 50;
    }
    if (m_fraction1 <= 0.0 || m_fraction1 > 1.0) {
        LogError(ERR_FILE_IO, "Invalid PADDS PerturbationValue. Defaulting to 0.2.");
        m_fraction1 = 0.2;
    }

    m_num_dec = m_pParamGroup->GetNumParams();
    m_num_objs = m_pObjFunc->CalcMultiObjFunc(NULL, -1);
    m_seed = GetRandomSeed();

    // Lower and upper bounds of the decision variables
    m_S_min.clear();
    m_S_max.clear();
    for (int i = 0; i < m_num_dec; i++) {
        m_S_min.push_back(m_pParamGroup->GetParamPtr(i)->GetLowerBoundTransformed());
        m_S_max.push_back(m_pParamGroup->GetParamPtr(i)->GetUpperBoundTransformed());
    }

    // Configure the multi-objective solution and initialize the workers
    SetMultiObjective(m_num_objs);
    std::vector<int> solveCounts(numberOfMpiProcesses, 0);
    ConfigureWorkers(solveCounts);

    // Write the information to the primary log file
    WriteSetup2(m_SetupName, m_ExecCmd, m_pObjFunc->GetObjFuncStr(), m_num_dec, m_pParamGroup->GetNumTiedParams());
    WriteStartingMetrics();
}

/**************************************************************************************************************************************************************
Optimize()

Search for the Pareto front using PADDS.
**************************************************************************************************************************************************************/
void PADDS::Optimize(void) {

    // Validate the setup, configure the workers and write the setup
    BeginAnalysis();

    // Number of initial solutions. User-specified initial solutions count towards the budget (per B. Tolson), and the random initial solutions are limited
    // to the remaining budget.
    int numInit = std::min((int)m_Init.size(), m_maxiter);
    int its = (int)(0.005 * m_maxiter);
    if (its < 5) its = 5;
    its = std::min(its, m_maxiter - numInit);

    // Create the initial solutions: the user-specified solutions followed by the random solutions
    std::vector<std::vector<double>> samples;
    for (int i = 0; i < numInit; i++) {
        samples.push_back(ApplyParameterGroup(m_Init[i]));
    }
    for (int i = 0; i < its; i++) {
        samples.push_back(RandomSolution());
    }

    // Solve the initial solutions together, then add them to the archive in order
    std::vector<std::vector<double>> objectives = Evaluate(samples);

    bool bRecorded = false;
    for (size_t i = 0; i < samples.size(); i++) {
        m_CurIter++;
        int result = UpdateArchive(samples[i], objectives[i]);
        bRecorded = false;
        if (result == ARCHIVE_NON_DOM) {
            WriteArchiveRecord(m_CurIter, (double)(m_maxiter - m_CurIter));
            bRecorded = true;
        }
    }

    // Number of evaluations remaining for the DDS search
    int iLeft = m_maxiter - its - numInit;

    // Calculating selection metric Z
    Calc_Z();

    // Main loop
    FILE* pPnFile = fopen("OstPADDSPn.txt", "w");
    fprintf(pPnFile, "EVAL  Pn\n");
    fclose(pPnFile);

    for (int i = 1; i < iLeft + 1; i++) {
        if (IsQuit() == true) { break; }

        // Select the solution to perturb
        Solution sbest = SelectParent();

        // DDS perturbation probability
        double Pn = 1 - log10((double)i) / log10((double)iLeft);

        pPnFile = fopen("OstPADDSPn.txt", "a");
        fprintf(pPnFile, "%04d  %E\n", i, Pn);
        fclose(pPnFile);

        // Perturb the selected solution and solve the candidate
        Solution stest = PerturbSolution(sbest, Pn);
        stest.F = Evaluate(std::vector<std::vector<double>>{ stest.X })[0];

        // Add the candidate to the archives
        bRecorded = ProcessCandidate(stest, sbest);
    }

    // Terminate the secondary workers and let them clean up.
    TerminateWorkers();

    // Write the final non-dominated archive, if it was not written after the last evaluation
    if (bRecorded == false) {
        WriteArchiveRecord(m_maxiter, 0.00);
    }

    // Write algorithm metrics
    WriteEndingMetrics();

} /* end Optimize() */

/**************************************************************************************************************************************************************
RandomSolution()

Generate a parameter set by sampling each decision variable uniformly within its (transformed) bounds
**************************************************************************************************************************************************************/
std::vector<double> PADDS::RandomSolution(void) {

    std::vector<double> sample;
    for (int j = 0; j < m_num_dec; j++) {
        sample.push_back(m_S_min[j] + (m_S_max[j] - m_S_min[j]) * UniformRandom());
    }
    return ApplyParameterGroup(sample);
}

/**************************************************************************************************************************************************************
SelectParent()

Select the solution to perturb. If the last candidate was dominated by its parent, a non-dominated solution is selected using the selection metric;
otherwise the most recently archived non-dominated solution is used.
**************************************************************************************************************************************************************/
PADDS::Solution PADDS::SelectParent(void) {

    if (m_dominance_flag == -1) {
        return SelectFrom();
    }
    return m_NonDom[0];
}

/**************************************************************************************************************************************************************
PerturbSolution()

Generate a candidate by perturbing each decision variable of the selected solution with probability Pn (at least one variable is perturbed)
**************************************************************************************************************************************************************/
PADDS::Solution PADDS::PerturbSolution(Solution& sbest, double Pn) {

    int dvn_count = 0;
    Solution stest = sbest;

    for (int j = 0; j < m_num_dec; j++) {
        if (UniformRandom() < Pn) {
            dvn_count += 1;
            stest.X[j] = neigh_value_continuous(sbest.X[j], m_S_min[j], m_S_max[j], m_fraction1);
        }
    }
    if (dvn_count == 0) {
        int dec_var = (int)(ceil(m_num_dec * UniformRandom() - 1.00));
        dec_var = std::max(0, std::min(dec_var, m_num_dec - 1));
        stest.X[dec_var] = neigh_value_continuous(sbest.X[dec_var], m_S_min[dec_var], m_S_max[dec_var], m_fraction1);
    }
    stest.X = ApplyParameterGroup(stest.X);

    return stest;
}

/**************************************************************************************************************************************************************
ProcessCandidate()

Add a solved candidate to the archives, given the solution it was generated from. Returns true if a new non-dominated solution was found and the run
record was written.
**************************************************************************************************************************************************************/
bool PADDS::ProcessCandidate(Solution& stest, Solution& sbest) {

    bool bRecorded = false;
    m_CurIter++;

    // Check to see if the selected solution dominates the candidate
    if (dominion_status(stest.F, sbest.F) == 2) {
        m_dominance_flag = -1;
        UpdateArchive(stest.X, stest.F);
    }
    else {
        // Check to see if the candidate duplicates a non-dominated solution. If so, discard it.
        if (IsDuplicate(stest.F)) {
            m_dominance_flag = 0;
            m_NumDuplicates++;
        }
        else {
            int result = UpdateArchive(stest.X, stest.F);
            if (result == ARCHIVE_NON_DOM) {
                WriteArchiveRecord(m_CurIter, (double)(m_maxiter - m_CurIter));
                bRecorded = true;
            }
        }

        if (m_dominance_flag != -1) {
            Calc_Z();
        }
    }

    return bRecorded;
}

/**************************************************************************************************************************************************************
Evaluate()

Solve a set of parameter sets and return the objectives of each
**************************************************************************************************************************************************************/
std::vector<std::vector<double>> PADDS::Evaluate(std::vector<std::vector<double>> samples) {

    std::vector<std::vector<double>> objectives(samples.size(), std::vector<double>(m_num_objs, INFINITY));
    if (samples.size() > 0) {
        ManageMultiObjectiveIterations(samples, 0, objectives);
    }
    return objectives;
}

/**************************************************************************************************************************************************************
ApplyParameterGroup()

Set a parameter set (transformed values) into the parameter group and read back the values. Integer parameters are rounded, so that the archived values
are the values that will be evaluated.
**************************************************************************************************************************************************************/
std::vector<double> PADDS::ApplyParameterGroup(std::vector<double> parameters) {

    SetParameterGroup(parameters);
    for (size_t entryParameter = 0; entryParameter < parameters.size(); entryParameter++) {
        parameters[entryParameter] = m_pParamGroup->GetParamPtr((int)entryParameter)->GetEstimatedValueTransformed();
    }
    return parameters;
}

/**************************************************************************************************************************************************************
SetParameterGroup()

Set a parameter set (transformed values) into the parameter group
**************************************************************************************************************************************************************/
void PADDS::SetParameterGroup(std::vector<double> parameters) {

    for (size_t entryParameter = 0; entryParameter < parameters.size(); entryParameter++) {
        m_pParamGroup->GetParamPtr((int)entryParameter)->SetEstimatedValueTransformed(parameters[entryParameter]);
    }
}

/**************************************************************************************************************************************************************
UpdateArchive()

Update the dominated and non-dominated archives with the latest solution. A new solution that is not dominated by any non-dominated solution is inserted at
the front of the non-dominated archive, and any non-dominated solutions that it (weakly) dominates are moved to the front of the dominated archive.
Otherwise the new solution is inserted at the front of the dominated archive.
**************************************************************************************************************************************************************/
int PADDS::UpdateArchive(std::vector<double> X, std::vector<double> F) {

    Solution newSolution;
    newSolution.X = X;
    newSolution.F = F;
    newSolution.Z = -999.999;

    // First entry is always non-dominated
    if (m_NonDom.empty() && m_Dom.empty()) {
        m_NonDom.push_back(newSolution);
        return ARCHIVE_NON_DOM;
    }

    // Assume the solution is non-dominated until we discover otherwise
    bool bMarkForInsertion = true;
    std::vector<Solution> remaining;

    // Compare against the current list of non-dominated solutions
    for (size_t entry = 0; entry < m_NonDom.size(); entry++) {
        Solution& current = m_NonDom[entry];

        // Does the new solution dominate the existing solution?
        bool bDominates = true;
        for (int i = 0; i < m_num_objs; i++) {
            if (current.F[i] < F[i]) {
                bDominates = false;
                break;
            }
        }

        // Is the new solution dominated by the existing solution?
        bool bIsDominated = false;
        if (bDominates == false) {
            bIsDominated = true;
            for (int i = 0; i < m_num_objs; i++) {
                if (F[i] < current.F[i]) {
                    bIsDominated = false;
                    break;
                }
            }
        }

        if (bDominates == true) {
            // Existing solution is dominated. Move it to the front of the dominated list.
            m_Dom.insert(m_Dom.begin(), current);
        }
        else {
            // Keep the existing solution
            remaining.push_back(current);

            // New solution is dominated. Make note so that it is not inserted into the list of non-dominated solutions.
            if (bIsDominated == true) {
                bMarkForInsertion = false;
            }
        }
    }
    m_NonDom = remaining;

    // Insert the new solution
    if (bMarkForInsertion == true) {
        m_NonDom.insert(m_NonDom.begin(), newSolution);
        return ARCHIVE_NON_DOM;
    }
    else {
        m_Dom.insert(m_Dom.begin(), newSolution);
        return ARCHIVE_DOM;
    }
}

/**************************************************************************************************************************************************************
IsDuplicate()

Check whether a non-dominated solution has exactly the same objectives
**************************************************************************************************************************************************************/
bool PADDS::IsDuplicate(std::vector<double>& F) {

    for (size_t entry = 0; entry < m_NonDom.size(); entry++) {
        if (m_NonDom[entry].F == F) {
            return true;
        }
    }
    return false;
}

/**************************************************************************************************************************************************************
dominion_status()

Determine the dominance relationship between two solutions. Returns 1 if F1 is no worse than F2 in every objective, 2 if F2 is no worse than F1 in every
objective (and F1 is worse in at least one), and 0 otherwise.
**************************************************************************************************************************************************************/
int PADDS::dominion_status(std::vector<double>& F1, std::vector<double>& F2) {

    bool bFirst = true;
    for (int i = 0; i < m_num_objs; i++) {
        if (F1[i] > F2[i]) {
            bFirst = false;
            break;
        }
    }
    if (bFirst) {
        return 1;
    }

    for (int i = 0; i < m_num_objs; i++) {
        if (F1[i] < F2[i]) {
            return 0;
        }
    }
    return 2;
}

/**************************************************************************************************************************************************************
SelectFrom()

Select a non-dominated solution with a probability proportional to its selection metric
**************************************************************************************************************************************************************/
PADDS::Solution PADDS::SelectFrom(void) {

    int archive_size = (int)m_NonDom.size();
    std::vector<double> z_cum(archive_size, 0.0);
    for (int i = 0; i < archive_size; i++) {
        if (i == 0) {
            z_cum[i] = m_NonDom[i].Z;
        }
        else {
            z_cum[i] = z_cum[i - 1] + m_NonDom[i].Z;
        }
    }

    double t = UniformRandom() * z_cum[archive_size - 1];

    int ii = 0;
    for (int i = 0; i < archive_size; i++) {
        if (z_cum[i] >= t) {
            ii = i;
            break;
        }
    }

    return m_NonDom[ii];
}

/**************************************************************************************************************************************************************
SortArchive()

Sort a list of indices into the non-dominated archive by the given objective (or by Z if whichObj == -1)
**************************************************************************************************************************************************************/
void PADDS::SortArchive(std::vector<int>& order, int whichObj) {

    int size = (int)order.size();
    for (int i = 0; i < size; i++) {
        for (int j = (i + 1); j < size; j++) {
            double F1, F2;
            if (whichObj == -1) {
                F1 = m_NonDom[order[i]].Z;
                F2 = m_NonDom[order[j]].Z;
            }
            else {
                F1 = m_NonDom[order[i]].F[whichObj];
                F2 = m_NonDom[order[j]].F[whichObj];
            }
            if (F2 < F1) {
                std::swap(order[i], order[j]);
            }
        }
    }
}

/**************************************************************************************************************************************************************
Calc_Z()

Compute the selection metric (Z) of each non-dominated solution
**************************************************************************************************************************************************************/
void PADDS::Calc_Z(void) {

    int archive_size = (int)m_NonDom.size();
    if (archive_size == 0) {
        return;
    }

    // Indices into the non-dominated archive, in archive order. These are sorted in place below.
    std::vector<int> pSorted(archive_size);
    for (int i = 0; i < archive_size; i++) {
        pSorted[i] = i;
    }

    switch (m_Select_metric) {
        case 0: { // Random
            for (int i = 0; i < archive_size; i++) {
                m_NonDom[i].Z = 1.0;
            }
            break;
        }
        case 1: { // Crowding distance
            for (int i = 0; i < archive_size; i++) {
                m_NonDom[i].Z = 0.0;
            }
            for (int obj = 0; obj < m_num_objs; obj++) {
                SortArchive(pSorted, obj);

                for (int i = 1; i < archive_size - 1; i++) {
                    double F0 = m_NonDom[pSorted[0]].F[obj];
                    double F1 = m_NonDom[pSorted[i - 1]].F[obj];
                    double F2 = m_NonDom[pSorted[i + 1]].F[obj];
                    double F4 = m_NonDom[pSorted[archive_size - 1]].F[obj];
                    m_NonDom[pSorted[i]].Z += fabs(F1 - F2) / fabs(F0 - F4);
                }

                if (archive_size > 1)
                    m_NonDom[pSorted[0]].Z = m_NonDom[pSorted[1]].Z;
                if (archive_size > 2)
                    m_NonDom[pSorted[archive_size - 1]].Z = m_NonDom[pSorted[archive_size - 2]].Z;
            }
            break;
        }
        case 2: { // Hypervolume contribution (Monte Carlo estimate)
            for (int i = 0; i < archive_size; i++) {
                m_NonDom[pSorted[i]].Z = 0;
            }

            // Define the boundaries
            std::vector<double> f_low_bound(m_num_objs), f_high_bound(m_num_objs);
            for (int i = 0; i < m_num_objs; i++) {
                SortArchive(pSorted, i);
                f_low_bound[i] = m_NonDom[pSorted[0]].F[i];
                f_high_bound[i] = m_NonDom[pSorted[archive_size - 1]].F[i];
            }

            // Prepare the Monte Carlo points
            int dots_num = 100;
            std::vector<std::vector<double>> mc_points(dots_num, std::vector<double>(m_num_objs));
            for (int i = 0; i < dots_num; i++) {
                for (int j = 0; j < m_num_objs; j++) {
                    mc_points[i][j] = f_low_bound[j] + (f_high_bound[j] - f_low_bound[j]) * UniformRandom();
                }
            }

            for (int i = 0; i < dots_num; i++) {
                int jj = archive_size;
                bool any_good = false;

                for (int j = 0; j < archive_size; j++) {
                    if (dominion_status(m_NonDom[pSorted[j]].F, mc_points[i]) == 1) {
                        jj = j;
                        any_good = true;
                        break;
                    }
                }

                if (any_good == true) {
                    for (int k = jj + 1; k < archive_size; k++) {
                        if (dominion_status(m_NonDom[pSorted[k]].F, mc_points[i]) == 1) {
                            // This point means nothing, go to the next point
                            any_good = false;
                            break;
                        }
                    }
                }

                // This point is only dominated by archive[jj]
                if (any_good == true) {
                    m_NonDom[pSorted[jj]].Z += 1.00;
                }
            }

            // Normalize Z
            double best_z = 0.00;
            for (int i = 0; i < archive_size; i++) {
                m_NonDom[pSorted[i]].Z = (m_NonDom[pSorted[i]].Z / (double)dots_num);
                if (m_NonDom[pSorted[i]].Z > best_z) {
                    best_z = m_NonDom[pSorted[i]].Z;
                }
            }
            for (int i = 0; i < archive_size; i++) {
                if (m_NonDom[pSorted[i]].Z == 0) {
                    m_NonDom[pSorted[i]].Z = 0.5 * best_z;
                }
            }
            break;
        }
        case 3: { // Hypervolume contribution (exact)
            int dataNumber = archive_size;
            int dimension = m_num_objs;

            std::vector<double> refPoint(dimension);
            for (int i = 0; i < dimension; i++) {
                SortArchive(pSorted, i);
                refPoint[i] = 1.00001 * (m_NonDom[pSorted[archive_size - 1]].F[i]);
            }

            // Hypervolume of the archive
            std::vector<std::vector<double>> pointsStore(dataNumber);
            std::vector<double*> pointsInitial(dataNumber);
            for (int i = 0; i < dataNumber; i++) {
                pointsStore[i] = m_NonDom[pSorted[i]].F;
                pointsInitial[i] = pointsStore[i].data();
            }
            double HyperVolume = HV(dataNumber, dimension, refPoint.data(), pointsInitial.data());

            // Contribution of each solution: the hypervolume lost when the solution is left out
            double best_z = 0;
            for (int i = 0; i < dataNumber; i++) {
                std::vector<double*> pointsInitial_sub;
                for (int j = 0; j < dataNumber; j++) {
                    if (j != i) {
                        pointsInitial_sub.push_back(pointsStore[j].data());
                    }
                }

                m_NonDom[pSorted[i]].Z = HyperVolume - HV(dataNumber - 1, dimension, refPoint.data(), pointsInitial_sub.data());

                if (m_NonDom[pSorted[i]].Z > best_z) {
                    best_z = m_NonDom[pSorted[i]].Z;
                }
            }

            // Taking care of the edges
            for (int i = 0; i < dimension; i++) {
                SortArchive(pSorted, i);
                m_NonDom[pSorted[0]].Z = best_z;
                m_NonDom[pSorted[archive_size - 1]].Z = best_z;
            }
            break;
        }
    }
}

/**************************************************************************************************************************************************************
neigh_value_continuous()

Perturb a decision variable using a normal distribution with a standard deviation of r times the variable range. Values outside of the bounds are
reflected or set to the bound with equal probability.
**************************************************************************************************************************************************************/
double PADDS::neigh_value_continuous(double s, double s_min, double s_max, double r) {

    double s_range = s_max - s_min;
    double snew = s + GaussRandom() * r * s_range;
    double P_Abs_or_Ref = UniformRandom();

    if (snew < s_min) {
        if (P_Abs_or_Ref <= 0.5) {
            snew = s_min + (s_min - snew);
        }
        else {
            snew = s_min;
        }
        if (snew > s_max) {
            snew = s_min;
        }
    }
    else if (snew > s_max) {
        if (P_Abs_or_Ref <= 0.5) {
            snew = s_max - (snew - s_max);
        }
        else {
            snew = s_max;
        }
        if (snew < s_min) {
            snew = s_max;
        }
    }
    return snew;
}

/**************************************************************************************************************************************************************
HV()

Hypervolume calculation.
**************************************************************************************************************************************************************/
double PADDS::HV(int data_n, int dim_n, double* ref, double** points) {

    int i, j;
    int dataNumber = data_n;

    // Set global
    m_dimension = dim_n;

    double** pointsInitial = new double* [dataNumber];
    for (int n = 0; n < dataNumber; n++) {
        pointsInitial[n] = new double[dim_n];
        for (int i = 0; i < dim_n; i++) {
            pointsInitial[n][i] = points[n][i];
        }
    }
    double* refPoint = new double[dim_n];
    for (i = 0; i < dim_n; i++) {
        refPoint[i] = ref[i];
    }

    // Initialize volume
    m_volume = 0.0;
    // Sqrt of dataNumber
    m_dSqrtDataNumber = sqrt((double)dataNumber);

    // Initialize region
    double* regionLow = new double[dim_n - 1];
    double* regionUp = new double[dim_n - 1];
    for (j = 0; j < dim_n - 1; j++) {
        // Determine minimal j coordinate
        double min = NEARLY_HUGE;
        for (i = 0; i < dataNumber; i++) {
            if (pointsInitial[i][j] < min) {
                min = pointsInitial[i][j];
            }
        }
        regionLow[j] = min;
        regionUp[j] = refPoint[j];
    }

    // Sort pointList according to d-th dimension
    SortPoints(pointsInitial, dataNumber, dim_n - 1);

    // Call stream initially
    stream(regionLow, regionUp, pointsInitial, dataNumber, 0, refPoint[m_dimension - 1]);

    for (int n = 0; n < dataNumber; n++) {
        delete[] pointsInitial[n];
    }
    delete[] pointsInitial;
    delete[] refPoint;
    delete[] regionLow;
    delete[] regionUp;

    return m_volume;
}

/**************************************************************************************************************************************************************
SortPoints()

Sort the array based on the which index.
**************************************************************************************************************************************************************/
void PADDS::SortPoints(double** X, int size, int which) {

    for (int i = 0; i < size; i++) {
        for (int j = (i + 1); j < size; j++) {
            if (X[j][which] < X[i][which]) {
                double* pTmp = X[i];
                X[i] = X[j];
                X[j] = pTmp;
            }
        }
    }
}

/**************************************************************************************************************************************************************
bool_vec_to_ulong()

Convert bitwise array of bools into equivalent integer: bit[0]*2^0 + bit[1]*2^1 + bit[2]*2^2 ...
**************************************************************************************************************************************************************/
int PADDS::bool_vec_to_ulong(bool* pB, int size) {

    int sum = 0;
    for (int i = 0; i < size; i++) {
        sum += (int)(pB[i] * pow((double)2, i));
    }
    return sum;
}

/**************************************************************************************************************************************************************
ulong_to_bool_vec()

Convert integer into equivalent bitwise array of bools.
**************************************************************************************************************************************************************/
void PADDS::ulong_to_bool_vec(int val, bool* pB, int size) {

    int mask = 1;
    for (int i = 0; i < size; i++) {
        pB[i] = ((mask & val) != 0);
        mask = mask << 1;
    }
}

/**************************************************************************************************************************************************************
covers()
**************************************************************************************************************************************************************/
bool PADDS::covers(double* cub, double* regLow) {

    for (int i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) {
        if (cub[i_hvc] > regLow[i_hvc]) {
            return false;
        }
    }
    return true;
}

/**************************************************************************************************************************************************************
partCovers()
**************************************************************************************************************************************************************/
bool PADDS::partCovers(double* cub, double* regUp) {

    for (int i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) {
        if (cub[i_hvc] >= regUp[i_hvc]) {
            return false;
        }
    }
    return true;
}

/**************************************************************************************************************************************************************
containsBoundary()
**************************************************************************************************************************************************************/
int PADDS::containsBoundary(double* cub, double* regLow, int split) {

    // Condition only checked for split>0
    if (regLow[split] >= cub[split]) {
        // Boundary in dimension split not contained in region, thus boundary is no candidate for the splitting line
        return -1;
    }
    else {
        for (int j_hvc = 0; j_hvc < split; j_hvc++) {
            if (regLow[j_hvc] < cub[j_hvc]) {
                // Boundary contained in region
                return 1;
            }
        }
    }
    // No boundary contained in region
    return 0;
}

/**************************************************************************************************************************************************************
getMeasure()
**************************************************************************************************************************************************************/
double PADDS::getMeasure(double* regLow, double* regUp) {

    double vol_hvc = 1.0;
    for (int i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) {
        vol_hvc *= (regUp[i_hvc] - regLow[i_hvc]);
    }
    return vol_hvc;
}

/**************************************************************************************************************************************************************
isPile()
**************************************************************************************************************************************************************/
int PADDS::isPile(double* cub, double* regLow, double* regUp) {

    int pile_hvc = m_dimension;

    // Check all dimensions of the node
    for (int k_hvc = 0; k_hvc < m_dimension - 1; k_hvc++) {
        // k-boundary of the node's region contained in the cuboid?
        if (cub[k_hvc] > regLow[k_hvc]) {
            if (pile_hvc != m_dimension) {
                // Second dimension occured that is not completely covered ==> cuboid is no pile
                return -1;
            }
            pile_hvc = k_hvc;
        }
    }

    // Region in only one dimension not completly covered ==> cuboid is a pile
    return pile_hvc;
}

/**************************************************************************************************************************************************************
computeTrellis()
**************************************************************************************************************************************************************/
double PADDS::computeTrellis(double* regLow, double* regUp, double* trellis) {

    int i_hvc, j_hvc;
    double vol_hvc = 0.0;
    double summand_hvc = 0.0;
    int numberSummands_hvc = 0;
    bool* bitvector_hvc = new bool[m_dimension];

    // Calculate number of summands
    bool* nSummands = new bool[m_dimension];
    for (i_hvc = 0; i_hvc < m_dimension; i_hvc++) {
        nSummands[i_hvc] = (i_hvc < m_dimension - 1);
    }
    numberSummands_hvc = bool_vec_to_ulong(nSummands, m_dimension);

    double* valueTrellis = new double[m_dimension - 1];
    double* valueRegion = new double[m_dimension - 1];
    for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) {
        valueTrellis[i_hvc] = trellis[i_hvc] - regUp[i_hvc];
    }
    for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) {
        valueRegion[i_hvc] = regUp[i_hvc] - regLow[i_hvc];
    }

    double* dTemp = new double[numberSummands_hvc / 2 + 1];

    // Sum
    for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) {
        // TODO Warning: dimension-1 <= 16 is assumed
        ulong_to_bool_vec(i_hvc, bitvector_hvc, m_dimension);

        // Construct summand. 0: take factor from region, 1: take factor from cuboid
        summand_hvc = 1.0;
        for (j_hvc = 0; j_hvc < m_dimension - 2; j_hvc++) {
            if (bitvector_hvc[j_hvc]) {
                summand_hvc *= valueTrellis[j_hvc];
            }
            else {
                summand_hvc *= valueRegion[j_hvc];
            }
        }
        summand_hvc *= valueRegion[m_dimension - 2];

        // Determine sign of summand
        vol_hvc -= summand_hvc;
        dTemp[i_hvc] = -summand_hvc;
    }

    ulong_to_bool_vec(i_hvc, bitvector_hvc, m_dimension);

    summand_hvc = 1.0;
    for (j_hvc = 0; j_hvc < m_dimension - 1; j_hvc++) {
        if (bitvector_hvc[j_hvc]) {
            summand_hvc *= valueTrellis[j_hvc];
        }
        else {
            summand_hvc *= valueRegion[j_hvc];
        }
    }
    vol_hvc -= summand_hvc;

    for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) {
        summand_hvc = dTemp[i_hvc];
        summand_hvc *= regUp[m_dimension - 2] - trellis[m_dimension - 2];
        summand_hvc /= valueRegion[m_dimension - 2];
        vol_hvc -= summand_hvc;
    }

    delete[] valueTrellis;
    delete[] valueRegion;
    delete[] dTemp;
    delete[] bitvector_hvc;
    delete[] nSummands;

    return vol_hvc;
}

/**************************************************************************************************************************************************************
getMedian()

Return median of the list of boundaries considered as a set
**************************************************************************************************************************************************************/
double PADDS::getMedian(double* bounds, int size) {

    // Do not filter duplicates
    if (size == 1) {
        return bounds[0];
    }
    else if (size == 2) {
        return bounds[1];
    }

    std::vector<double> pSort(bounds, bounds + size);
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (pSort[j] < pSort[i]) {
                std::swap(pSort[i], pSort[j]);
            }
        }
    }

    return pSort[size / 2];
}

/**************************************************************************************************************************************************************
stream()

Recursive calculation of hypervolume
**************************************************************************************************************************************************************/
void PADDS::stream(double* regionLow, double* regionUp, double** points, int npoints, int split, double cover) {

    //--- init --------------------------------------------------------------//
    double coverOld_hvc = cover;
    unsigned int coverIndex_hvc = 0;
    int c_hvc;

    //--- cover -------------------------------------------------------------//

    // Identify first covering cuboid
    double dMeasure = getMeasure(regionLow, regionUp);
    while (cover == coverOld_hvc && coverIndex_hvc < (unsigned int)npoints) {
        if (covers(points[coverIndex_hvc], regionLow)) {
            // New cover value
            cover = points[coverIndex_hvc][m_dimension - 1];
            m_volume += dMeasure * (coverOld_hvc - cover);
        }
        else coverIndex_hvc++;
    }

    // coverIndex shall be the index of the first point in points which is ignored in the remaining process. It may occur that that some points in front
    // of coverIndex have the same d-th coordinate as the point at coverIndex. These points must be discarded and therefore the following loop checks for
    // these points and reduces coverIndex if necessary.
    for (c_hvc = coverIndex_hvc; c_hvc > 0; c_hvc--) {
        if (points[c_hvc - 1][m_dimension - 1] == cover) {
            coverIndex_hvc--;
        }
    }

    // Abort if points is empty
    if (coverIndex_hvc == 0) {
        return;
    }
    // Note: in the remainder points is only considered to index coverIndex

    //--- allPiles  ---------------------------------------------------------//

    bool allPiles = true;
    unsigned int iii;

    std::vector<int> piles_hvc(coverIndex_hvc);
    for (iii = 0; iii < coverIndex_hvc; iii++) {
        piles_hvc[iii] = isPile(points[iii], regionLow, regionUp);
        if (piles_hvc[iii] == -1) {
            allPiles = false;
            break;
        }
    }

    // trellis[i] contains the values of the minimal i-coordinate of the i-piles. If there is no i-pile the default value is the upper bound of the region.
    // The 1-dimensional KMP of the i-piles is: reg[1][i] - trellis[i]
    if (allPiles) { // sweep

        // Initialize trellis with region's upper bound
        std::vector<double> trellis(regionUp, regionUp + (m_dimension - 1));

        double current = 0.0;
        double next = 0.0;
        iii = 0;
        do { // while(next != coverNew)
            current = points[iii][m_dimension - 1];
            do { // while(next == current)
                if (points[iii][piles_hvc[iii]] < trellis[piles_hvc[iii]]) {
                    trellis[piles_hvc[iii]] = points[iii][piles_hvc[iii]];
                }
                iii++; // index of next point
                if (iii < coverIndex_hvc) {
                    next = points[iii][m_dimension - 1];
                }
                else {
                    next = cover;
                }

            } while (next == current);
            m_volume += computeTrellis(regionLow, regionUp, trellis.data()) * (next - current);
        } while (next != cover);
    }
    //--- split -------------------------------------------------------------//
    // Inner node of partition tree
    else {
        double bound = -1.0;
        int boundaries_size = 0;
        int noBoundaries_size = 0;
        std::vector<double> boundaries(coverIndex_hvc);
        std::vector<double> noBoundaries(coverIndex_hvc);

        do {
            for (iii = 0; iii < coverIndex_hvc; iii++) {
                int contained = containsBoundary(points[iii], regionLow, split);
                if (contained == 1) {
                    boundaries[boundaries_size] = points[iii][split];
                    boundaries_size++;
                }
                else if (contained == 0) {
                    noBoundaries[noBoundaries_size] = points[iii][split];
                    noBoundaries_size++;
                }
            }

            if (boundaries_size > 0) {
                bound = getMedian(boundaries.data(), boundaries_size);
            }
            else if (noBoundaries_size > m_dSqrtDataNumber) {
                bound = getMedian(noBoundaries.data(), noBoundaries_size);
            }
            else {
                split++;
            }
        } while (bound == -1.0);

        double dLast;
        int pointsChild_size = 0;
        std::vector<double*> pointsChild(coverIndex_hvc);

        // Left child: reduce maxPoint
        dLast = regionUp[split];
        regionUp[split] = bound;
        for (iii = 0; iii < coverIndex_hvc; iii++) {
            if (partCovers(points[iii], regionUp)) {
                pointsChild[pointsChild_size] = points[iii];
                pointsChild_size++;
            }
        }
        if (pointsChild_size > 0) {
            stream(regionLow, regionUp, pointsChild.data(), pointsChild_size, split, cover);
        }

        // Right child: increase minPoint
        pointsChild_size = 0;

        regionUp[split] = dLast;
        dLast = regionLow[split];
        regionLow[split] = bound;
        for (iii = 0; iii < coverIndex_hvc; iii++) {
            if (partCovers(points[iii], regionUp)) {
                pointsChild[pointsChild_size] = points[iii];
                pointsChild_size++;
            }
        }
        if (pointsChild_size > 0) {
            stream(regionLow, regionUp, pointsChild.data(), pointsChild_size, split, cover);
        }
        regionLow[split] = dLast;
    }
}

/**************************************************************************************************************************************************************
GetSelectionMetricString()

Get the description of the selection metric
**************************************************************************************************************************************************************/
const char* PADDS::GetSelectionMetricString(void) {

    if (m_Select_metric == 0) return "random";
    else if (m_Select_metric == 1) return "crowding distance";
    else if (m_Select_metric == 2) return "estimated hypervolume contribution";
    else if (m_Select_metric == 3) return "exact hypervolume contribution";
    return "unknown";
}

/**************************************************************************************************************************************************************
WriteHeader()

Write the column names: one objective per observation group, followed by the parameters
**************************************************************************************************************************************************************/
void PADDS::WriteHeader(FILE* pFile, const char* pAfter) {

    for (int entryObjective = 0; entryObjective < m_num_objs; entryObjective++) {
        fprintf(pFile, "obj.%-10s  ", m_pObsGroup->GetGroup(entryObjective));
    }
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "%s\n", pAfter);
}

/**************************************************************************************************************************************************************
WriteSolution()

Write the objectives and parameters of a solution (without a line ending)
**************************************************************************************************************************************************************/
void PADDS::WriteSolution(FILE* pFile, Solution& solution) {

    for (int entryObjective = 0; entryObjective < m_num_objs; entryObjective++) {
        fprintf(pFile, "%E  ", solution.F[entryObjective]);
    }
    SetParameterGroup(solution.X);
    m_pParamGroup->Write(pFile, WRITE_SCI);
}

/**************************************************************************************************************************************************************
WriteArchiveRecord()

Write the non-dominated archive to the run record (OstOutput0.txt), the console and OstNonDomSolutions0.txt
**************************************************************************************************************************************************************/
void PADDS::WriteArchiveRecord(int iter, double remaining) {

    FILE* pFiles[2];
    pFiles[0] = fopen("OstOutput0.txt", "a");
    pFiles[1] = stdout;

    for (FILE* pFile : pFiles) {
        fprintf(pFile, "\n");
        for (size_t entry = 0; entry < m_NonDom.size(); entry++) {
            fprintf(pFile, "%-4d  ", iter);
            WriteSolution(pFile, m_NonDom[entry]);
            fprintf(pFile, "%E\n", remaining);
        }
    }
    fclose(pFiles[0]);

    // The current non-dominated solutions
    FILE* pFile = fopen("OstNonDomSolutions0.txt", "w");
    fprintf(pFile, "gen   ");
    WriteHeader(pFile, "trials remaining");
    for (size_t entry = 0; entry < m_NonDom.size(); entry++) {
        fprintf(pFile, "%-4d  ", iter);
        WriteSolution(pFile, m_NonDom[entry]);
        fprintf(pFile, "%E\n", remaining);
    }
    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteArchives()

Write the final non-dominated and dominated archives
**************************************************************************************************************************************************************/
void PADDS::WriteArchives(FILE* pFile) {

    fprintf(pFile, "\nNon-Dominated Solutions\n");
    WriteHeader(pFile, "");
    for (size_t entry = 0; entry < m_NonDom.size(); entry++) {
        WriteSolution(pFile, m_NonDom[entry]);
        fprintf(pFile, "\n");
    }

    fprintf(pFile, "\nDominated Solutions\n");
    WriteHeader(pFile, "");
    for (size_t entry = 0; entry < m_Dom.size(); entry++) {
        WriteSolution(pFile, m_Dom[entry]);
        fprintf(pFile, "\n");
    }

    fprintf(pFile, "\nNumber of Non-Dominated Solutions : %d\n", (int)m_NonDom.size());
    fprintf(pFile, "\nNumber of Dominated Solutions     : %d\n", (int)m_Dom.size());
}

/**************************************************************************************************************************************************************
WriteStartingMetrics()

Write out the PADDS setup and the header of the run record
**************************************************************************************************************************************************************/
void PADDS::WriteStartingMetrics(void) {

    FILE* pFile = fopen("OstOutput0.txt", "a");

    fprintf(pFile, "Max Iterations          : %d\n", m_maxiter);
    fprintf(pFile, "Perturbation Value      : %lf\n", m_fraction1);
    fprintf(pFile, "Selection Metric        : %s\n", GetSelectionMetricString());
    fprintf(pFile, "Num Objectives          : %d\n", m_num_objs);
    fprintf(pFile, "Initial Param. Sets     : %d\n", (int)m_Init.size());
    if (m_BatchSize > 0) {
        fprintf(pFile, "Batch Size              : %d\n", m_BatchSize);
    }

    fprintf(pFile, "gen   ");
    WriteHeader(pFile, "trials remaining");

    fclose(pFile);
}

/**************************************************************************************************************************************************************
WriteEndingMetrics()

Write out the final archives and the algorithm metrics
**************************************************************************************************************************************************************/
void PADDS::WriteEndingMetrics(void) {

    // Final archives, to the run record and the console
    FILE* pFile = fopen("OstOutput0.txt", "a");
    WriteArchives(pFile);
    WriteArchives(stdout);

    // Algorithm metrics
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : %s\n", m_AlgorithmName.c_str());
    fprintf(pFile, "Max Iterations          : %d\n", m_maxiter);
    fprintf(pFile, "Actual Iterations       : %d\n", m_CurIter);
    fprintf(pFile, "Num Decision Variables  : %d\n", m_num_dec);
    fprintf(pFile, "Num Objectives          : %d\n", m_num_objs);
    fprintf(pFile, "Random Seed             : %u\n", m_seed);
    fprintf(pFile, "Perturbation Value      : %lf\n", m_fraction1);
    fprintf(pFile, "Non-Dominated Solutions : %d\n", (int)m_NonDom.size());
    fprintf(pFile, "Dominated Solutions     : %d\n", (int)m_Dom.size());
    fprintf(pFile, "Duplicates Discarded    : %d\n", m_NumDuplicates);
    fprintf(pFile, "Selection Metric        : %s\n", GetSelectionMetricString());
    if (m_BatchSize > 0) {
        fprintf(pFile, "Batch Size              : %d\n", m_BatchSize);
    }
    fprintf(pFile, "Total Evals             : %d\n", m_NumSolves);

    fclose(pFile);
}

/**************************************************************************************************************************************************************
PADDS_Program()

Search for the Pareto front using PADDS.
**************************************************************************************************************************************************************/
void PADDS_Program(int argC, StringType argV[]) {

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Determine behavior based on the rank
    if (rank == 0) {
        // This is the primary worker. Setup the algorithm and call the processing routine
        PADDS* pPADDS = new PADDS();

        // Begin the primary worker operations. This transfers information to the secondary workers and manages the solution process
        pPADDS->Optimize();
    }
    else {
        // This is a secondary worker. Create the secondary worker, telling it to configure from MPI, and start the work
        ModelWorker worker = ModelWorker(true);
        worker.WorkMPI();
    }
} /* end PADDS_Program() */
