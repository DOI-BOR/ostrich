/******************************************************************************
File     : StatUtility.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

StatUtility contains C routines that assist in statistical calculations.

Version History
06-19-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   replaced magic # with NEARLY_ZERO
01-01-07    lsm   Added Durbin-Watson and runs tests.
******************************************************************************/
#ifndef STAT_UTILITY_H
#define STAT_UTILITY_H

#include "MyHeaderInc.h"

/* what type of central tendency to use in statistical calculations */
#define CENTRAL_TEND_MEAN   (0)
#define CENTRAL_TEND_MEDIAN (1)
#define CENTRAL_TEND_PCTILE (2)

double GammaLn(double val);
double CalcMean(Ironclad1DArray  v, int size);
double CalcMedian(Unmoveable1DArray v, int size);
double CalcStdDev(Ironclad1DArray v, int size, int ctType);
double CalcSkewness(Ironclad1DArray v, int size);
double CalcKurtosis(Ironclad1DArray v, int size);
double GetCritValNormPPCC(int n);

double FdistCDF(int u, int v, double xLwr, double xUpr);
double FdistPDF(int u, int v, double x);
double FdistInvCDF(int u, int v, double p);

double StudentCDF(int dof, double x);
double StudentPDF(int dof, double x);
double StudentInvCDF(int dof, double p);

double StdNormCDF(double x);
double StdNormPDF(double x);
double StdNormInvCDF(double p);

bool RunsTest(double * residuals, int num, int * nPos, int * nNeg, int * nRuns, 
              int * clwr, int * cupr);
double RunsTestCDF(int nRuns, int nPos, int nNeg);
int InvRunsTestCDF(double a, int m, int n);

void AutorunFunctionTest(double * residuals, int num, double * r1, double * var, 
                         double * vpx, double * med, int * nSur, int * nDef, 
                         int * n1, double * clwr, double * cupr);

double nCr(int n, int r);

void STATS_TestFdist(void);
void STATS_TestStudentDist(void);
void STATS_TestStdNormDist(void);
void STATS_TestAutocorrelation(void);

#endif /* STAT_UTILITY_H */

