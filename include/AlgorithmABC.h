/******************************************************************************
File     : AlgorithmABC.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

AlgorithmABC is a base class which defines the common interface to be used by 
the set of optimization/calibration algorithms that make up the Ostrich program.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added MAX_COUNT 
08-13-04    lsm   converted to ABC
******************************************************************************/
#ifndef ALGORITHM_ABC_H
#define ALGORITHM_ABC_H

#include "MyHeaderInc.h"

#define MAX_COUNT (3) //max. # iterations to tolerate without dec. obj. func.

/******************************************************************************
class AlgorithmABC
   Algorithms utilize the Optimize() subroutine to perform optimization on a 
   given Objective Function, or use Calibrate() to perform Least-Squares
   regression on a given Model.
******************************************************************************/
class AlgorithmABC
{
   public:
	  virtual ~AlgorithmABC(void) { DBG_PRINT("AlgorithmABC::DTOR"); }
      virtual void Destroy(void)=0;
      virtual void Optimize(void)=0;
      virtual void Calibrate(void)=0;      
      virtual void WriteMetrics(FILE * pFile)=0;
      virtual void WarmStart(void)=0;
      virtual int  GetCurrentIteration(void)=0;
}; /* end class Algorithm */

#endif /* ALGORITHM_ABC_H */

