/******************************************************************************
File     : Exception.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

Defines the various kinds of errors that can occur along with a convenient 
C-style interface for reporting and recovering from such errors.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-24-04    lsm   added MemCheck() function.
07-08-04    lsm   added Jacobian errror, added GetParameterByName()
08-17-04    lsm   Now limiting the number of error messages stored in the log. 
                  User abort error was added, added GetTiedParameterByName()
12-10-04    lsm   Error file is now deleted prior to program execution
10-19-05    lsm   Added support for BGA and GRID programs
01-01-07    lsm   Added five more error codes:
                     ERR_BAD_WGHT --> invalid observation weight
                     ERR_INS_PARM --> insensitive parameter
                     ERR_INS_OBS  --> insensitive observation
                     ERR_CONTINUE --> the error message is a continutation of a previous msg
07-13-07    lsm   Added SuperMUSE error code (ERR_SMUSE) and SuperMUSE cleanup.
******************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Model.h"
#include "AlgorithmABC.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "TelescopingBounds.h"
#include "StatsClass.h"

#include "Exception.h"
#include "SuperMuseUtility.h"
#include "Utility.h"
#include "IsoParse.h"

/*
Definition of an error message structure, containing: (1) an error code, 
(2) an error string describing the error, and (3) a pointer to the next 
error message.
*/
typedef struct ERROR_MSG_STRUCT
{
   ErrorCodeType             errCode;
   char                      errMsg[DEF_STR_SZ];
   struct ERROR_MSG_STRUCT * pNxt;
}ErrorMsg;

/* maximum number of errors that can be logged (to prevent memory overflow) */
#define MAX_ERRORS (1000)

/*
A mapping between error enums and human readable strings.
*/
IroncladString ErrorMap[NUM_ERRORS] = 
{
   "NO ERROR",
   "BAD ARGUMENTS",
   "FILE I/O ERROR",
   "MODEL EXECUTION ERROR",
   "ARRAY OUT OF BOUNDS",
   "PARAMETER MISMATCH",
   "SINGULAR MATRIX",
   "GRID SIZE IS TOO LARGE",
   "INITIAL SA TEMPERATURE",
   "PARAMETER BOUNDS",
   "COULDN'T BOUND MINIMUM",
   "UNKNOWN BOUND CONDITION",
   "COULDN'T PARSE INPUT",
   "MALLOC/NEW FAILED",
   "JACOBIAN ERROR",
   "USER ABORT",
   "BINARY CODED GA",
   "OBSERVATION WEIGHTS",
   "INSENSITIVE PARAMETER",
   "INSENSITIVE OBSERVATION",
   " ",
   "SUPERMUSE",
   "OVERFLOW (DIV-BY-ZERO?)",
   "NULL POINTER",
   "ALGORITHM STALLED",
   "FILE CLEANUP",
   "NON-UNIQUE PARAMETER NAME",
   "FIXED FORMAT PARAMETERS",
   "DEGREES OF FREEDOM"
};

/******************************************************************************
Exception Variables

These variables store information about any errors that occur in the 
program.
******************************************************************************/
ErrorMsg m_ErrList; //the list of errors
int gNumErrors = 0; //number of errors (includes ones not logged, due to overflow) */
char m_ErrorFile[DEF_STR_SZ];
int  m_Id;  //MPI processor id
//constructor and destructor counters
long int m_CTORS = 0;
long int m_DTORS = 0;
//registered pointers to model and algorithms objects, for exiting gracefully
ModelABC * m_pModelReg = NULL;
AlgorithmABC * m_pAlgReg = NULL;
StatsClass * m_pStatsReg = NULL;
char m_pIterResPrefix[1000];

//user-aborts program by creating this file
const char * gStopFile = "OstQuit.txt";

//behavioral threshold (used by PreserveModelOutput)
double gObjFuncThreshold = HUGE_VAL;

//trial number (used by PreserveModelOutput)
int gTrialNumber = 0;

void SetObjFuncThreshold(double threshold)
{
   gObjFuncThreshold = threshold;
}

double GetObjFuncThreshold(void)
{
   return(gObjFuncThreshold);
}

void SetTrialNumber(int trial_number)
{
   gTrialNumber = trial_number;
}

int GetTrialNumber(void)
{
   return(gTrialNumber);
}

//clean up error list
void DestroyErrorList(ErrorMsg * pMsg);

/******************************************************************************
RunModel()

Execute the model.
******************************************************************************/
double RunModel(void)
{
   double f = HUGE_VAL;
   if(m_pModelReg != NULL)
   {
      f = m_pModelReg->Execute();     
   }
   return f;
}/* end RunModel() */

/******************************************************************************
SetIterationResidualsPrefix()

Set a prefix to use for the iteration residuals file name. This is used to 
idetify the trial (iPrefix) and algorithm (sPrefix) for the residuals in a 
multi-start or hybrid algorithm.
******************************************************************************/
void SetIterationResidualsPrefix(char * sPrefix, int iPrefix)
{
   char tmpStr[1000];

   if (sPrefix != NULL)
   {
      strcpy(m_pIterResPrefix, sPrefix);
   }
   else
   {
      sprintf(tmpStr, "_T%03d", iPrefix);
      strcat(m_pIterResPrefix, tmpStr);
   }   
} /* SetIterationResidualsPrefix */

/******************************************************************************
WriteIterationResiduals()

Save iterration residuals to file.
******************************************************************************/
void WriteIterationResiduals(void)
{
   int step;

   if (m_pStatsReg == NULL) return;
   if (m_pModelReg == NULL) return;
   if (m_pAlgReg == NULL) return;

   step = m_pAlgReg->GetCurrentIteration();

   m_pStatsReg->WriteResiduals(step, m_pIterResPrefix);
}/* end WriteIterationResiduals() */

/******************************************************************************
SaveModel()

Save the model files.
******************************************************************************/
void SaveModel(int id)
{
   if(m_pModelReg != NULL)
   {
      m_pModelReg->SaveBest(id);
   }
}/* end SaveModel() */

/******************************************************************************
RegisterModelPtr() and RegisterAlgPtr()

Registrs model and algorithm pointers. When program fails, the registered 
classes are destroyed so that the abnormal program termination does not cause 
a memory leak.
******************************************************************************/
void RegisterModelPtr(ModelABC * pModel) { m_pModelReg = pModel; strcpy(m_pIterResPrefix, "");}
void RegisterAlgPtr(AlgorithmABC * pAlg) { m_pAlgReg = pAlg;}

/******************************************************************************
RegisterStatsPtr()

Registers a pointer to the statistics member class of an algorithm. This allows
for logging of residuals at each step of the algorithm.
******************************************************************************/
void RegisterStatsPtr(StatsClass * pStats){ m_pStatsReg = pStats;}

/******************************************************************************
TelescopicCorrection()

Correct candidate parameter value using telescoping strategy.
******************************************************************************/
double TelescopicCorrection(double xmin, double xmax, double xbest, double a, 
                           double xnew)
{
   if(m_pModelReg == NULL) return xnew;

   switch(m_pModelReg->GetTelescopingStrategy())
   {
      case (TSCOPE_NONE) : return xnew;
      case (TSCOPE_PVEX) : return (telescope_parameter(xmin, xmax, xbest, a, xnew, &fpvx));
      case (TSCOPE_CVEX) : return (telescope_parameter(xmin, xmax, xbest, a, xnew, &fvex));
      case (TSCOPE_LINR) : return (telescope_parameter(xmin, xmax, xbest, a, xnew, &flin));
      case (TSCOPE_CAVE) : return (telescope_parameter(xmin, xmax, xbest, a, xnew, &fcve));
      case (TSCOPE_DCVE) : return (telescope_parameter(xmin, xmax, xbest, a, xnew, &fdcv));
      default : return xnew;
   }/* end switch() */
}/* end TelescopicCorrection() */

/******************************************************************************
GetParameterByName() and GetTiedParameterByName()

Retrieves pointers to Paramters and TiedParameters which have the name 
specified in the input argument. Returns NULL if no such name is found.
******************************************************************************/
ParameterABC * GetParameterByName(IroncladString pName)
{
   if(m_pModelReg == NULL){ return NULL;}
   if(m_pModelReg->GetParamGroupPtr() == NULL){ return NULL;}
   return(m_pModelReg->GetParamGroupPtr()->GetParamPtr(pName));
}/* end GetParameterByName() */

/******************************************************************************
GetTiedParameterByName()

Retrieves ptr to tied parameter with given name.
******************************************************************************/
TiedParamABC * GetTiedParameterByName(IroncladString pName)
{
   if(m_pModelReg == NULL){ return NULL;}
   if(m_pModelReg->GetParamGroupPtr() == NULL){ return NULL;}
   return(m_pModelReg->GetParamGroupPtr()->GetTiedParamPtr(pName));
}/* end GetTiedParameterByName() */

/******************************************************************************
GetConstraintByName()

Retrieves ptr to contraint with given name.
******************************************************************************/
ConstraintABC * GetConstraintByName(IroncladString pName)
{
   if(m_pModelReg == NULL){ return NULL;}
   if(m_pModelReg->GetObjFuncPtr() == NULL){ return NULL;}
   return(m_pModelReg->GetObjFuncPtr()->GetConstraintPtr(pName));
}/* end GetConstraintByName() */

/******************************************************************************
GetParameterName()

Retrieves name of Paramter at given index, or NULL if not initialized.
******************************************************************************/
IroncladString GetParameterName(int idx)
{
   if(m_pModelReg == NULL){ return NULL;}
   if(m_pModelReg->GetParamGroupPtr() == NULL){ return NULL;}
   return(m_pModelReg->GetParamGroupPtr()->GetParamPtr(idx)->GetName());
}/* end GetParameterName() */

/******************************************************************************
GetParameterValStr()

Retrieves a string representation of the value of the Paramter at the given 
index.
******************************************************************************/
IroncladString GetParameterValStr(int idx, double val)
{
   static char str[DEF_STR_SZ];
   double old;

   if(m_pModelReg == NULL){ return NULL;}
   if(m_pModelReg->GetParamGroupPtr() == NULL){ return NULL;}

   old = m_pModelReg->GetParamGroupPtr()->GetParamPtr(idx)->GetEstVal();
   m_pModelReg->GetParamGroupPtr()->GetParamPtr(idx)->SetEstVal(val);
   m_pModelReg->GetParamGroupPtr()->GetParamPtr(idx)->GetValAsStr(str);
   m_pModelReg->GetParamGroupPtr()->GetParamPtr(idx)->SetEstVal(old);

   return(str);
}/* end GetParameterValStr() */

/******************************************************************************
GetNumDigitsOfPrecision()

Retrieves the number of digits of precision to use when writing outputs or 
model inputs.
******************************************************************************/
int GetNumDigitsOfPrecision(void)
{
   if(m_pModelReg == NULL){ return 6;}
   return(m_pModelReg->GetNumDigitsOfPrecision());
}/* end GetNumDigitsOfPrecision() */

/******************************************************************************
GetFilePairs()

Retrieves the list of file pairs for a given model.
******************************************************************************/
FilePair * GetFilePairs(void)
{
   if(m_pModelReg == NULL){ return NULL;}
   return(((Model *)m_pModelReg)->GetFilePairs());
}/* end GetFilePairs() */

/******************************************************************************
IncCtorCount() and IncDtorCount()

Increments Constructor and Destructor counters. When the program finishes, the
number of DTOR calls should be equal to the number of CTOR calls. Otherwise,
there may be a memory leak.
******************************************************************************/
void IncCtorCount(void) {m_CTORS++;}
void IncDtorCount(void) {m_DTORS++;}

/******************************************************************************
InitErrors()

Initializes variables to default values.
******************************************************************************/
void InitErrors(void)
{
   MPI_Comm_rank(MPI_COMM_WORLD, &m_Id);
   sprintf(m_ErrorFile, "OstErrors%d.txt", m_Id);
   
   strcpy(m_ErrList.errMsg, "No messages");
   m_ErrList.errCode = ERR_NO_ERROR;
   m_ErrList.pNxt = NULL;

   remove(gStopFile);
   remove(m_ErrorFile);
}/* end InitErrors() */

/******************************************************************************
SetErrorFile()

Overrides the default output filename.
******************************************************************************/
void SetErrorFile(IroncladString filename)
{
   strcpy(m_ErrorFile, filename);
}/* end SetErrorFile() */

/******************************************************************************
ReportErrors()

Writes error information to standard output and also to file. Once written,
the logged messages are freed from RAM.
******************************************************************************/
void ReportErrors(void)
{
   UnchangeableString errStr;   
   FILE * pFile;
   ErrorMsg * pErr;

   pFile = fopen(m_ErrorFile, "w");
   fprintf(stdout,"Ostrich Error Report for Processor %d \n", m_Id);
   fprintf(pFile, "Ostrich Error Report for Processor %d \n", m_Id);

   fprintf(stdout,"A total of %d errors and/or warnings were reported\n", gNumErrors);
   fprintf(pFile, "A total of %d errors and/or warnings were reported\n", gNumErrors);

   pErr = &m_ErrList;
   while(pErr != NULL)
   {
      errStr = ErrorMap[(int)pErr->errCode];      
      fprintf(stdout, "%-26s : ", errStr);   
      fprintf(stdout, "%s \n", pErr->errMsg);
      
      fprintf(pFile, "%-26s : ", errStr);
      fprintf(pFile, "%s \n", pErr->errMsg);

      pErr = pErr->pNxt;
   } /* end while() */

   if(gNumErrors > MAX_ERRORS)
   {
      fprintf(stdout,"Warning: The number of errors/warnings (%d) exceeded the max size of the error list (%d)\n", gNumErrors, MAX_ERRORS);
      fprintf(pFile, "Warning: The number of errors/warnings (%d) exceeded the max size of the error list (%d)\n", gNumErrors, MAX_ERRORS);
      
      fprintf(stdout,"Only the first %d errors/warnings were logged in the error file.\n", MAX_ERRORS);
      fprintf(pFile, "Only the first %d errors/warnings were logged in the error file.\n", MAX_ERRORS);
   }

   DestroyErrorList(m_ErrList.pNxt);
   m_ErrList.pNxt = NULL;
   gNumErrors = 0;

   fclose(pFile);
}/* end ReportErrors() */

/******************************************************************************
DestroyErrorList()

Free up error list.
******************************************************************************/
void DestroyErrorList(ErrorMsg * pMsg)
{
   if(pMsg != NULL)
   {
      DestroyErrorList(pMsg->pNxt);
      delete pMsg;
   }
}/* end DestroyErrorList() */

/******************************************************************************
GetErrorCode()

Retrieves the last error code.
******************************************************************************/
ErrorCodeType GetErrorCode(void)
{
   ErrorMsg * pErr;
   pErr = &m_ErrList;

   while(pErr->pNxt != NULL)
   {
      pErr = pErr->pNxt;
   }/* end while() */

   return (pErr->errCode);
}/* end GetErrorCode() */

/******************************************************************************
LogError()

Adds the error code and message to the list.
******************************************************************************/
void LogError(ErrorCodeType err, IroncladString msg)
{
   ErrorMsg * pErr;
   pErr = &m_ErrList;

   gNumErrors++;
   //limit number of errors that can be logged
   if(gNumErrors <= MAX_ERRORS)
   {
      if(pErr->errCode != ERR_NO_ERROR)
      {
         while(pErr->pNxt != NULL)
         {
            pErr = pErr->pNxt;
         }/* end while() */

         NEW_PRINT("ErrorMsg", 1);
         pErr->pNxt = new ErrorMsg;
         MEM_CHECK(pErr->pNxt);

         pErr = pErr->pNxt;
      }/* end if() */

      strncpy(pErr->errMsg, msg, DEF_STR_SZ);
      // log messages are constrained to be no longer than DEF_STR_SZ
      if(strlen(msg) > DEF_STR_SZ) 
      {
         pErr->errMsg[DEF_STR_SZ-1] = NULLSTR;
      }
      pErr->errCode = err;
      pErr->pNxt = NULL;
   }/* end if() */
}/* end LogError() */

/******************************************************************************
FileOpenFailure()

Reports a file open error and closes down the program.
******************************************************************************/
void FileOpenFailure(IroncladString routine, IroncladString file)
{
   char msg[DEF_STR_SZ];
   
   //set error and exit program
	sprintf(msg, "%s(): couldn't open |%s|\n", routine, file);
   LogError(ERR_FILE_IO, msg);   
   ExitProgram(1);
}/* end FileOpenFailure() */

/******************************************************************************
EndOfFileFailure()

Reports that an input file ended unexpectedly and closes down the program.
******************************************************************************/
void EndOfFileFailure(IroncladString routine, IroncladString file)
{
   char msg[DEF_STR_SZ];

   //set error and exit program
	sprintf(msg, "%s(): %s file input ended unexpectedly ", routine, file);
   LogError(ERR_FILE_IO, msg);   
   ExitProgram(1);
}/* end EndOfFileFailure() */

/******************************************************************************
MissingTokenFailure()

Reports that an input file is missing a token.
******************************************************************************/
void MissingTokenFailure(IroncladString token, IroncladString file)
{
   char msg[DEF_STR_SZ];

   //set error and exit program
	sprintf(msg, "Missing token %s in file %s", token, file);
   LogError(ERR_FILE_IO, msg);   
   ExitProgram(1);
} /* end MissingTokenFailure() */

/******************************************************************************
ExitProgram()

Quits the program, but gracefully.
******************************************************************************/
void ExitProgram(int code)
{   
   ReportErrors();

   if(code != 0)
   {
      delete m_pAlgReg;
      delete m_pModelReg;

      MPI_Abort(MPI_COMM_WORLD, 1);   
   }  

   DestroySuperMUSE();

   fprintf(stdout, "num CTORS: %ld \n", m_CTORS);
   fprintf(stdout, "num DTORS: %ld \n", m_DTORS);

   //free up copy array, used in matrix inversion
   MatInv(NULL, NULL, 0);

   //free up Ostrich data line
   InitDataLine(NULL);

   //free up IsoFit data line
   ISO_GetFileSize(NULL);

   //delete temporary input file
   #ifndef ISOFIT_BUILD
      if(strcmp(GetOstFileName(), GetInFileName()) != 0)
      {
      	remove(GetInFileName());
      }
   #endif

   exit(0);
}/* end ExitProgram() */

/******************************************************************************
MemCheck()

Checks memory.
******************************************************************************/
void MemCheck(void * pMem, int line, const char * file)
{
   char Msg[DEF_STR_SZ];

   if(pMem == NULL)
   { 
      sprintf(Msg, "Memory allocation error on line %d of file %s!",line,file); 
      LogError(ERR_MALLOC, Msg);
      ExitProgram(1);
   }
}/* end MemCheck() */

/******************************************************************************
IsQuit()

Checks to see if the user has requested early termination of the program.
******************************************************************************/
bool IsQuit(void)
{
   static bool isLogged = false;
   FILE * pFile;
   pFile = fopen(gStopFile, "r");
   if(pFile != NULL)
   { 
      fclose(pFile);
      if(isLogged == false)
      {
         LogError(ERR_ABORT, "stop file detected, aborting program");
         isLogged = true;
      }
      return true;
   }

   return false;
}/* end IsQuit() */

