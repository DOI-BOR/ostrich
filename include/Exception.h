/******************************************************************************
File     : Exception.h
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

This file defines a c-style interface for tracking and reporting any errors 
that cause the program to halt prematurely.

At the beginning of main(), InitErrors() must be called.

Prior to returning from main(), the contents of the error message and error 
code should be ouput to the user via a call to ReportErrors().

If a critical error occurs in the program, a call to LogError() should be 
made prior to exiting.

ReportErrors() will write errors to OstErrors.txt unless SetErrorFile() is 
called after InitErrors() is called.

FileOpenFailure() can be called whenever a file fails to open. This routine
will report the offending routine and file along with the current working 
directory, and then terminate the program.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-26-03    lsm   created version history field and updated comments.
11-19-03    lsm   Added WRITE_TX_BNR definition.
07-08-04    lsm   added Jacobian errror, added WRITE_OPT, doubled LINE_SIZE,
                  added GetParameterByName()
08-17-04    lsm   Now limiting the number of error messages stored in the log. 
                  User abort error was added, added GetTiedParameterByName()
10-19-05    lsm   Added support for BGA and GRID programs
01-01-07    lsm   Added five more error codes:
                     ERR_BAD_WGHT --> invalid observation weight
                     ERR_INS_PARM --> insensitive parameter
                     ERR_INS_OBS  --> insensitive observation
                     ERR_CONTINUE --> the error message is a continutation of a previous msg
07-13-07    lsm   Added SuperMUSE error code (ERR_SMUSE).
******************************************************************************/
#ifndef MY_EXCEPTION_H
#define MY_EXCEPTION_H

#include "MyHeaderInc.h"

//forward class declarations
class ModelABC;
class AlgorithmABC;
class StatsClass;
class ParameterABC;
class TiedParamABC;
class ConstraintABC;
class FilePair;

#define ERR_VALUE (-999.999)

#define WRITE_SCI (1) //write values in scientific format (%E)
#define WRITE_DEC (2) //write values in decimal format (%.6lf)
#define WRITE_BNR (3) //write names/banner of values (%s)
//write names/banner of values (%s) also add information about transformations
#define WRITE_TX_BNR (4) 
#define WRITE_DBG (5) //write all variables
#define WRITE_OPT (6) //write using the 'optimal' format

/* 
An enum is defined for any error that can cause the program to halt
prematurely.
*/
typedef enum ERROR_CODE_TYPE
{
   ERR_NO_ERROR = 0,
   ERR_BAD_ARGS = 1,
   ERR_FILE_IO  = 2,
   ERR_MODL_EXE = 3,  /* execution of the model (split) */
   ERR_ARR_BNDS = 4,  /* array out of bounds */
   ERR_MISMATCH = 5,  /* parameter mismatch */
   ERR_SING_MAT = 6,  /* singular matrix */
   ERR_GRD_SIZE = 7,  /* grid size too large */
   ERR_SA_TEMP  = 8,  /* initial simulated annealing temperature */
   ERR_PRM_BNDS = 9,  /* parameter outside of bounds */
   ERR_BND_MIN  = 10, /* could not bound min */
   ERR_BND_UNK  = 11, /* unknown min bounding condition */  
   ERR_IN_PARSE = 12, /* failed to parse a line of input */
   ERR_MALLOC   = 13, /* couldn't allocate memory */
   ERR_JACOBIAN = 14, /* Jacobian insensitivity */
   ERR_ABORT    = 15, /* user abort */
   ERR_BGA      = 16, /* error in binary-coded GA */
   ERR_BAD_WGHT = 17, /* invalid observation weight */
   ERR_INS_PARM = 18, /* insensitive parameter */
   ERR_INS_OBS  = 19, /* insensitive observation */
   ERR_CONTINUE = 20, /* the error message is a continutation of a previous msg */
   ERR_SMUSE    = 21, /* error related to SuperMUSE parallel cluster */
   ERR_OVERFLOW = 22, /* overflow condition (typically caused by divide-by-zero) */
   ERR_NULL_PTR = 23, /* NULL pointer detected */
   ERR_STALL    = 24, /* algorithm stall detected */
   ERR_CLEANUP  = 25,
   ERR_PRM_NEST = 26,  /* parameter is substring of another parameter */
   ERR_FIXD_FMT = 27,  /* error in conversion to fixed format */
   ERR_BAD_DOF  = 28   /* insufficient degrees of freedom */
}ErrorCodeType;

/* NUM_ERRORS is used to size ErrorMap defined in Exception.cpp */
#define NUM_ERRORS (29) 

/* 
These c-style routines allow access to the unique instance of the 
error structure defined in Exception.cpp.
*/
extern "C" 
{
   double TelescopicCorrection(double xmin, double xmax, double xbest, double a, double xnew);
   bool IsQuit(void);
   void InitErrors(void);
   void ReportErrors(void);   
   ErrorCodeType GetErrorCode(void);
   void SetErrorFile(IroncladString filename);
   void LogError(ErrorCodeType err, IroncladString msg);
   void FileOpenFailure(IroncladString routine, IroncladString file);
   void EndOfFileFailure(IroncladString routine, IroncladString file);
   void MissingTokenFailure(IroncladString token, IroncladString file);
   void ExitProgram(int code);

   void RegisterModelPtr(ModelABC * pModel);
   void RegisterAlgPtr(AlgorithmABC * pAlg);
   void RegisterStatsPtr(StatsClass * pStats);

   int GetNumDigitsOfPrecision(void);
   ParameterABC * GetParameterByName(IroncladString pName);
   TiedParamABC * GetTiedParameterByName(IroncladString pName);
   ConstraintABC * GetConstraintByName(IroncladString pName);
   IroncladString GetParameterName(int idx);
   IroncladString GetParameterValStr(int idx, double val);
   FilePair * GetFilePairs(void);
   #define MEM_CHECK(a) MemCheck((void *)(a), __LINE__, __FILE__);

   void MemCheck(void * pMem, int line, const char * file);
   double RunModel(void);
   void SaveModel(int id);
   void WriteIterationResiduals(void);
   void SetIterationResidualsPrefix(char * sPrefix, int iPrefix);

   void SetObjFuncThreshold(double threshold);
   double GetObjFuncThreshold(void);

   void SetTrialNumber(int trial_number);
   int GetTrialNumber(void);
}

#endif /* MY_EXCEPTION_H */

