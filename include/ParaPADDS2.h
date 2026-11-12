/**************************************************************************************************************************************************************

Parallel PADDS (Pareto Archived Dynamically Dimensioned Search). This is the parallel version of the PADDS multi-objective algorithm, and shares the PADDS
archive, selection metric and hypervolume functions.

The pre-refactor implementation was asynchronous: each worker was sent a new candidate as soon as it returned a result. The Algorithm base class solves
sets of alternatives using a blocking batch call, so this version is synchronous. Candidates are generated in batches (one per available worker) from the
current archive, solved concurrently, and then added to the archive in the order they were generated. Each candidate is compared against the solution it
was generated from.

**************************************************************************************************************************************************************/

#ifndef PARA_PADDS_ALGORITHM_H
#define PARA_PADDS_ALGORITHM_H

// Include custom classes
#include "PADDS2.h"


class ParaPADDS : public PADDS {
    public:
        // Define constructors and destructors
        ParaPADDS();
        ~ParaPADDS(void) { DBG_PRINT("ParaPADDS::DTOR"); }

        // Solution functions
        void Optimize(void);                                                // Function to search for the Pareto front using parallel PADDS

    private:
        void WriteFrontObjectives(void);                                    // Function to write the objectives of the non-dominated solutions (c_results.txt)
};

extern "C" {
void PARA_PADDS_Program(int argC, StringType argV[]);
}

#endif /* PARA_PADDS_ALGORITHM_H */
