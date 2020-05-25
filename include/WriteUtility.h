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

// forward decs
class ModelABC;
class AlgorithmABC;

extern "C" {
void WritePreciseNumber(FILE * pOut, double x);
void WriteDisclaimer(FILE * pFile);
void WriteSetup(ModelABC * pModel, const char * algStr);
void WriteSetupToFile(FILE * pFile, ModelABC * pModel, const char * algStr);
void WriteSetupNoDisclaimer(ModelABC * pModel, const char * algStr);

void WriteBanner(ModelABC * pModel, const char * pBef, const char * pAft);
void WriteBannerToFile(FILE * pFile, ModelABC * pModel, const char * pBef, const char * pAft);

void WriteMultiObjRecord(ModelABC * pModel, int iter, ArchiveStruct * pArch, double dx);
void WriteMultiObjRecordToFile(FILE * pFile, ModelABC * pModel, int iter, ArchiveStruct * pArch, double dx);

void WriteRecord(ModelABC * pModel, int iter, double fx, double dx);
void WriteRecordToFile(FILE * pFile, ModelABC * pModel, int iter, double fx, double dx);
void WriteStatus(StatusStruct * pStatus);

void WriteMultiObjOptimal(ModelABC * pModel, ArchiveStruct * pNonDom, ArchiveStruct * pDom);
void WriteMultiObjOptimalToFile(FILE * pFile, ModelABC * pModel, ArchiveStruct * pNonDom, ArchiveStruct * pDom);

void WriteOptimal(ModelABC * pModel, double fx);
void WriteOptimalToFile(FILE * pFile, ModelABC * pModel, double fx);

void WriteAlgMetrics(AlgorithmABC * pAlg);
void WriteMelt(int count, int max, char c);

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

void Write1dSearch(int count, int max);
void WriteInnerEval(int count, int max, char c);

void WriteGrid(GridStruct * pGrid, int size);
}
#endif /* WRITE_UTILITY_H */

