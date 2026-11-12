/******************************************************************************
File      : SCEUA2.h
Author    : L. Shawn Matott - Converted from Origincal SCEUA Fortran code
Copyright : 2009, L. Shawn Matott

The SCE-UA method is a general purpose global optimization
program.  It was originally developed by Dr. Qingyun Duan as part
of his doctoral dissertation work at the Department of Hydrology
and Water Resources, University of Arizona, Tucson, AZ 85721, USA. 
The dissertation is entitled "A Global Optimization Strategy for
Efficient and Effective Calibration of Hydrologic Models".  The
program has since been modified to make it easier for use on
problems of users' interests.  The algorithm has been described
in detail in an article entitled "Effective and Efficient Global
Optimization for Conceptual Rainfall-Runoff Models", Water Resources
Research, Vol 28(4), pp.1015-1031, 1992; and in an article entitled
"A Shuffled Complex Evolution Approach for Effective and Efficient
Global Minimization" by Q. Duan, V.K. Gupta and S. Sorooshian,
Journal of Optimization Theory and its Applications, Vol 76(3), 
pp 501-521, 1993.  A paper entitled "Optimal Use of the SCE-UA Global
Optimization Method for Calibrating Watershed Models", by Q. Duan,
S. Sorooshian, & V.K. Gupta, Journal of Hydrology, Vol.158, 265-284,
1994, discussed how to use the SCE-UA Method in an efficient and 
effective manner.

Input Summary for the SCEUA algorithm (adapted from original Fortran-based
description):
==========================================================================
variable   type     description
MAXN       integer  Maximum number of trials allowed before
                    optimization is terminated.  The purpose of
                    MAXN is to stop an optimization search before
                    too much computer time is expended.  MAXN
                    should be set large enough so that optimization
                    is generally completed before MAXN trials are
                    performed. Recommended value is 10,000 (increase or
                    decrease as necessary).
---> this parameter is called m_Budget within Ostrich

KSTOP      integer  Number of shuffling loops in which the 
                    criterion must improve by the specified
                    percentage or else optimization will be
                    terminated. Recommended value: 5.
---> this parameter is called m_Kstop within Ostrich

PCENTO     double   Percentage by which the criterion value must
                    change in the specified number of shuffling 
                    loops or else optimization is terminated
                    (Use decimal equivalent: Percentage/100).
                    Recommended value: 0.01.
---> this parameter is called m_Pcento within Ostrich

NGS        integer  Number of complexes used for optimization
                    search.  Minimum value is 1.
                    Recommended value is between 2 and 20 depending
                    on the number of parameters to be optimized and
                    on the degree of difficulty of the problem.
---> this parameter is called m_NumComplexes within Ostrich

ISEED      integer  Random seed used in optimization search.  Enter
                    any integer number.  Default value (=1969) is
                    assumed if this field is left blank.
                    Recommended value: any large integer.
---> this parameter is called m_Seed within Ostrich

IDEFLT     integer  Flag for setting the control variables of the
                    SCE-UA algorithm.  Enter false or leave the field
                    blank for default setting.  Enter true for user
                    specified setting.
                    If option true is chosen, user must specify alg. 
                    parameters.
---> this parameter is called m_bUserConfig within Ostrich

NPG        integer  Number of points in each complex.  NPG should
                    be greater than or equal to 2.  The default
                    value is equal to (2 * number of optimized
                    parameters + 1).
---> this parameter is called m_PtsPerComplex within Ostrich

NPS        integer  Number of points in each sub-complex.  NPS
                    should be greater than or equal to 2 and less
                    than NPG.  The default value is equal to 
                    (number of optimized parameters + 1).
---> this parameter is called m_PtsPerSubComplex within Ostrich

NSPL       integer  Number of evolution steps taken by each complex
                    before next shuffling.  Default value is equal
                    to NPG.
---> this parameter is called m_NumEvoSteps within Ostrich

MINGS      integer  Minimum number of complexes required for
                    optimization search, if the number of complexes
                    is allowed to reduce as the optimization search
                    proceeds.  The default value is equal to NGS.
---> this parameter is called m_MinComplexes within Ostrich

INIFLG     bool     Flag on whether to include an initial point in
                    the starting population.  Enter true if the initial 
                    point is to be included.  The default value is equal to false.
---> this parameter is called m_bUseInitPt within Ostrich

IPRINT    integer   Print-out control flag.  Enter '0' to print out
                    the best estimate of the global optimum at the
                    end of each shuffling loop.  Enter '1' to print
                    out every point in the entire sample population
                    at the end of each shuffling loop.  The default
                    value is equal to 0. Enter 2 to ignore this variable
                    and use conventional Ostrich output.
---> this parameter is called m_OutputMode within Ostrich

PARAMS     double * Initial estimates of the parameters to be optimized.
---> this parameter is called m_pParams within Ostrich

LOWER      double * Lower bounds of the parameters to be optimized.
---> this parameter is called m_pLower within Ostrich

UPPER      double * Upper bounds of the parameters to be optimized.
---> this parameter is called m_pUpper within Ostrich

Version History
10-31-09    lsm   Created
******************************************************************************/
#ifndef SCEUA_H
#define SCEUA_H

// Include C classes
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// Include C++ classes
#include <vector>
#include <iostream>

// Include custom classes
#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"


class SCEUA : public Algorithm {
   public:
        // Define constructors and destructors
        SCEUA();
        ~SCEUA(void){ DBG_PRINT("SCEUA::DTOR"); Destroy(); }
        void Destroy(void);

        // Warm restart is not yet supported for the SCEUA algorithm
        void WarmStart(void);

        // Solution functions
        void Optimize(void);                                                // Function to optimize a solution using SCEUA
        void Calibrate(void);                                               // Function to solve least squares fitting using SCEUA
        int  GetCurrentIteration(void) { return m_CurIter; }

        // Output functions
        void WriteStartingMetrics(void);                                    // Function to write the starting SCEUA metrics
        void WriteEndingMetrics(void);                                      // Function to write the ending SCEUA metrics

   private:
        // Functions converted from the original FORTRAN implementation
        void InitFromFile(IroncladString pFileName);
        void scemain(void);
        void scein(double * a, double * bl, double * bu, int nopt, int *maxn,
                   int *kstop, double *pcento, int *iseed, int *ngs, int *npg,
                   int *nps, int *nspl, int *mings, int *iniflg, int *iprint);
        void sceua(double * a, double * bl, double * bu, int nopt, int maxn,
                   int kstop, double pcento, int iseed, int ngs, int npg,
                   int nps, int nspl, int mings, int iniflg, int iprint);
        void cce(int nopt, int nps, double ** s, double * sf, double * bl,
                 double * bu, double * xnstd, int * icall, double maxn,
                 int * iseed);
        void getpnt(int nopt, int idist, int * iseed, double * x, double * bl,
                    double * bu, double * std, double * xi);
        void sort(int n, int m, double ** rb, double * ra);
        void parstt(int npt, int nopt, double ** x, double * xnstd,
                    double * bound, double * gnrng, int * ipcnvg);
        void sort(int n, int * ra);
        void comp(int n, int npt, int ngs1, int ngs2, int npg,
                  double ** a, double * af, double ** b, double * bf);
        void chkcst(int nopt, double * snew, double * bl,
                    double * bu, int * ibound);

        // Model evaluation functions
        void ApplyParameterCorrections(double * x);                         // Applies parameter corrections to a point and stores the corrected values back into it
        void EvaluateAlternatives(std::vector<std::vector<double>> alternatives,
                                  std::vector<double>& objectives);         // Solves a set of points across the workers and updates the running best
        double EvaluateAlternative(double * x, int nopt);                   // Solves a single point and updates the running best

        // Configuration variables
        double m_Best = INFINITY;                                           // Best objective of the current population
        int m_Budget; //MAXN
        int m_Kstop; //KSTOP
        double m_Pcento; //PCENTO
        double m_Peps; //peps
        int m_NumComplexes; //NGS
        int m_Seed; //ISEED
        int m_UserConfig; //IDEFLT
        int m_PtsPerComplex; //NPG
        int m_PtsPerSubComplex; //NPS
        int m_NumEvoSteps; //NSPL
        int m_MinComplexes; //MINGS
        int m_UseInitPt; //INIFLG
        int m_OutputMode; //IPRINT
        int m_np; //number of parameters
        double * m_pParams = NULL; //PARAMS
        double * m_pLower = NULL; //LOWER
        double * m_pUpper = NULL; //UPPER
        bool m_bUseInitPt = false;

        // Solution tracking variables
        int m_CurIter = 0;                                                  // Counter for the number of records written to the log file
        double m_PctComplete = 0.00;                                        // Percentage of the budget that has been used
}; /* end class SCEUA */

extern "C" {
void SCEUA_Program(int argC, StringType argV[]);
}

#endif /* SCEUA_H */
