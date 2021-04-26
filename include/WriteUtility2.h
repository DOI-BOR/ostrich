/******************************************************************************
File      : WriteUtility.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

This file provides a unifying interface for the various algorithms to write output
both to file and to stdout.

Version History
03-08-04    lsm   created
08-11-04    lsm   upped version, added date and time of build, added PATO and
                  algorithm metrics support
01-28-05    lsm   Added support for tracking algorithm status (used in 
                  grid-computing). Added WriteGrid()
01-01-07    lsm   upped version, added support for ModelABC
******************************************************************************/
#ifndef WRITE_UTILITY_H
#define WRITE_UTILITY_H

#include "MyHeaderInc.h"
#include "Algorithm.h"
#include "ParameterGroup.h"

// forward decs
class Algorithm;

extern "C" {
void WritePreciseNumber2(FILE * pOut, double x);
void WriteDisclaimer2(FILE * pFile);
void WriteSetup2(Algorithm *algo, const char * algStr);
void WriteSetupToFile2(FILE * pFile, Algorithm *algo, const char * algStr);
void WriteSetupNoDisclaimer2(Algorithm *algo, const char * algStr);

void WriteBanner2(Algorithm *algo, const char * pBef, const char * pAft);
void WriteBannerToFile2(FILE * pFile, Algorithm *algo, const char * pBef, const char * pAft);

void WriteMultiObjRecord2(Algorithm *algo, int iter, ArchiveStruct * pArch, double dx);
void WriteMultiObjRecordToFile2(FILE * pFile, Algorithm *algo, int iter, ArchiveStruct * pArch, double dx);

void WriteRecord2(Algorithm *algo, int iter, double fx, double dx);
void WriteRecordToFile2(FILE * pFile, Algorithm *algo, int iter, double fx, double dx);

void WriteMultiObjOptimal2(Algorithm *algo, ArchiveStruct * pNonDom, ArchiveStruct * pDom);
void WriteMultiObjOptimalToFile2(FILE * pFile, Algorithm *algo, ArchiveStruct * pNonDom, ArchiveStruct * pDom);

void WriteOptimal2(Algorithm *algo, double fx);
void WriteOptimalToFile2(FILE * pFile, Algorithm *algo, double fx);
void WriteOptimalToFileWithGroup2(FILE* pFile, ParameterGroup *paramGroup, double fx);

//void WriteAlgMetrics(AlgorithmABC * pAlg);
void WriteMelt2(int count, int max, char c);

#define WRITE_BRENT  (0)
#define WRITE_GSECT (-1)
#define WRITE_SWTCH (-2)
#define WRITE_ENDED (-3)
#define WRITE_GA    (-4)
#define WRITE_PSO   (-5)
#define WRITE_SA    (-6)
#define WRITE_LEV   (-7)
#define WRITE_GRID  (-8)
#define WRITE_BIS   (-9)
#define WRITE_SMP   (-10)
#define WRITE_DDS   (-11)
#define WRITE_LHS   (-12)
#define WRITE_USR   (-13)
#define WRITE_JAC   (-14)
#define WRITE_SCE   (-15)
#define WRITE_GLUE  (-16)

void Write1dSearch2(int count, int max);
void WriteInnerEval2(int count, int max, char c);

void WriteGrid2(GridStruct * pGrid, int size);
}
#endif /* WRITE_UTILITY_H */

