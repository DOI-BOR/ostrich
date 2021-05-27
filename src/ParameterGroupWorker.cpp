/******************************************************************************
File      : ParameterGroup.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates a group of parameters. The optimization routines will attempt to
find the values of the parameters that minimize the objective function.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-19-03    lsm   Modified to support multiple stages of unit conversion.
07-05-04    lsm   added integer and combinatorial parameter support
07-08-04    lsm   added tied parameter support
07-09-04    lsm   added CheckTemplateFiles() to check for parameters that are 
                  not in any template files
12-01-04    lsm   Added support for Geometry parameters. Also added support for
                  random initialization of parameters.
10-19-05    lsm   Added CheckBounds(). Replaced rand() with MyRand()
01-01-07    lsm   Added ExcludeParam() subroutine to support the "hold" 
                  parameters functionality. Added two new tied parameter types 
                  to support ratios: TiedParamSimpleRatio and TiedParamComplexRatio.
07-18-07    lsm   Added support for SuperMUSE
******************************************************************************/
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "ParameterGroupWorker.h"

#include "FilePair.h"
#include "FilePipe.h"
#include "ConstraintABC.h"
#include "DatabaseABC.h"
#include "VertexList.h"

#include "FortranSupportUtilities.h"
#include "Utility.h"
#include "Exception.h"

/******************************************************************************
GetNumParams()

Retrieves the number of parameters.
******************************************************************************/
int ParameterGroupWorker::GetNumParams(void) {
  return m_NumParams;
} 

/******************************************************************************
ReadParams()

Stuffs an array with the current parameter values. Array must have been 
previously allocated.
******************************************************************************/
void ParameterGroupWorker::ReadParams(double * p) {
   for(int j = 0; j < m_NumParams; j++) {
      p[j] = m_pList[j]->GetValue();
   }
}

/******************************************************************************
GetParamPtr()

Retrieves a pointer to the ith parameter.
******************************************************************************/
ParameterWorkerABC* ParameterGroupWorker::GetParamPtr(int i)
{
  return m_pList[i];
} /* end GetParamPtr() */

/******************************************************************************
GetParamPtr()

Retrieves a pointer to the parameter with matching name.
******************************************************************************/
ParameterWorkerABC* ParameterGroupWorker::GetParamPtr(IroncladString name)
{
   int j;

   //determine indices by examing parameter names
   for(j = 0; j < m_NumParams; j++)
   {
      if(strcmp(m_pList[j]->GetName(), name) == 0){ return m_pList[j];}
   }/* end for() */   
   return NULL;
} /* end GetParamPtr() */


/******************************************************************************
Destroy()

Frees the memory used by the objects of all the parameters contained in it
******************************************************************************/
void ParameterGroupWorker::Destroy(void)
{
   int j;
   for(j = 0; j < m_NumParams; j++)
   {
      delete m_pList[j];
   }
   delete [] m_pList;

   for(j = 0; j < m_NumExcl; j++)
   {
      delete m_pExcl[j];
   }
   delete [] m_pExcl;

   for(j = 0; j < m_NumParams; j++)
   {
      delete [] m_ParamNameList[j];
   }
   delete [] m_ParamNameList;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR

Initializes parameter group from user-specified input file.
******************************************************************************/
ParameterGroupWorker::ParameterGroupWorker() {
   m_pList = NULL;
   m_pExcl = NULL;
   m_ParamNameList = NULL;
   m_NumParams = 0;
   m_NumExcl = 0;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SubIntoFile()

Substitutes the estimated value of the parameter into the model input file.
******************************************************************************/
void ParameterGroupWorker::SubIntoFile(FilePipe * pPipe)
{ 
   int i, size;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   char * pRep;
   ParameterWorkerABC* pParam;


   //Adjustable parameters
   for(i = 0; i < m_NumParams; i++)
   {
      pParam = m_pList[i];
      strcpy(find, pParam->GetName());
      pParam->GetValAsStr(replace);
      pPipe->FindAndReplace(find,replace);
   } /* end for() */

   //Excluded parameters
   for(i = 0; i < m_NumExcl; i++)
   {
      pParam = m_pExcl[i];
      strcpy(find, pParam->GetName());
      pParam->GetValAsStr(replace);
      pPipe->FindAndReplace(find,replace);
   } /* end for() */

   pPipe->StringToFile();
} /* end SubIntoFile() */

/******************************************************************************
WriteDatabaseParameter()

Loop over database entries until the desired parameter is written.
******************************************************************************/
void ParameterGroupWorker::WriteDatabaseParameter(DatabaseABC * pDbase, char * find, char * replace)
{
   char ErrorMsg[DEF_STR_SZ];
   DatabaseABC * pCur = NULL;
   bool bFound = false;
   for(pCur = pDbase; pCur != NULL; pCur = pCur->GetNext())
   {
      if(pCur->WriteParameter(find,replace) == true) bFound = true; /* parameter could be in multiple databases */
      //if(pCur->WriteParameter(find,replace) == true)
      //   break;
   }
   if(bFound == false)
   {
      sprintf(ErrorMsg, "Parameter |%s| not found in list of database entries!", find);
      LogError(ERR_MISMATCH, ErrorMsg);
   }
}/* end WriteDatabaseParameter() */

/******************************************************************************
SubIntoDbase()

Substitutes the estimated value of the parameter into the model input database.
******************************************************************************/
void ParameterGroupWorker::SubIntoDbase(DatabaseABC * pDbase)
{ 
   int i, size;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   char * pRep;
   ParameterWorkerABC* pParam;

   //Adjustable parameters
   for(i = 0; i < m_NumParams; i++)
   {
      pParam = m_pList[i];
      strcpy(find, pParam->GetName());
      pParam->GetValAsStr(replace);
      WriteDatabaseParameter(pDbase, find, replace);
   } /* end for() */

   //Excluded parameters
   for(i = 0; i < m_NumExcl; i++)
   {
      pParam = m_pExcl[i];
      strcpy(find, pParam->GetName());
      pParam->GetValAsStr(replace);
      WriteDatabaseParameter(pDbase, find, replace);
   } /* end for() */

} 

/******************************************************************************
WriteSuperMuseArgs()

This routine is similar in purpose to SubIntoFile(), except that file substution
information (i.e. find/replace pairs) are written to a SuperMUSE arguments file.
Final substitution into the model template files will be performed by a SuperMUSE
client-side tasker script/batch file.
******************************************************************************/
void ParameterGroupWorker::WriteSuperMuseArgs(FILE * pFile)
{
   int i, size;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   char * pRep;
   ParameterWorkerABC* pParam;

   //Adjustable parameters
   for(i = 0; i < m_NumParams; i++)
   {
      pParam = m_pList[i];
      strcpy(find, pParam->GetName());
      pParam->GetValAsStr(replace);
      fprintf(pFile, "%s %s ", find, replace);
   } /* end for() */

   //Excluded parameters
   for(i = 0; i < m_NumExcl; i++)
   {
      pParam = m_pExcl[i];
      strcpy(find, pParam->GetName());
      pParam->GetValAsStr(replace);
      fprintf(pFile, "%s %s ", find, replace);
   } /* end for() */

   fprintf(pFile, " \n");
}/* end WriteSuperMuseArgs() */


/******************************************************************************
GetNextEmptyParamIdx()

Finds the next unassigned parameter in the array.
******************************************************************************/
int ParameterGroupWorker::GetNextEmptyParamIdx(void)
{
   for(int i = 0; i < m_NumParams; i++){ if(m_pList[i] == NULL) { return i;}}

   LogError(ERR_ARR_BNDS, "GetNextEmptyParamIdx() : array is filled!");   
   return 0;
}/* end GetNextEmptyParamIdx() */

/******************************************************************************
ExtractInitialValue()

Reads the parameter value from a model input file. Returns true if successful.
******************************************************************************/
bool ParameterGroupWorker::ExtractValue(char * name, bool bFixFmt, double * pVal)
{
   bool bFound;
   double val;
   //extract parameters from model input file
   FilePipe * pPipe;
   FilePair * pCur = GetFilePairs();
   while(pCur != NULL)
   {
      pPipe = pCur->GetPipe();
      bFound = ExtractParameter(name, pPipe->GetTemplateFileName(), 
                                pPipe->GetModelInputFileName(), bFixFmt, &val, m_ParamNameList, m_NumParams);
      if(bFound)
      {
         *pVal = val;
         return true;
      }
      pCur = pCur->GetNext();
   } /* end while() */

   return false;
}/* end ExtractInitialValue() */

/******************************************************************************
SetGroupValues()

Allows ability to set values into a parameter group
******************************************************************************/
void ParameterGroupWorker::SetGroupValues(ParameterWorkerABC** m_pListInput, ParameterWorkerABC** m_pExclInput, char** m_ParamNameListInput,
                                          int numberOfParameters,  int numberOfExcluded) {

    // Set the group values from the input arrays
    m_pList = m_pListInput;
    m_pExcl = m_pExclInput;
    m_ParamNameList = m_ParamNameListInput;

    // Set the size values from the input arrays
    m_NumParams = numberOfParameters;
    m_NumExcl = numberOfExcluded;
}

/******************************************************************************
Write()

Writes formatted output to the pFile argument.
******************************************************************************/
void ParameterGroupWorker::Write(FILE * pFile, int type)
{

   for(int i = 0; i < m_NumParams; i++) {      
      m_pList[i]->Write(pFile, type);
   }

} /* end writeToFile() */

/******************************************************************************
CheckTemplateFiles()

Checks to see if every parameter is included in at least one template file.
Parameters not found in any template file will trigger a warning message but 
will not halt the program.
******************************************************************************/
void ParameterGroupWorker::CheckTemplateFiles(FilePair * pList)
{
   FilePair * pCur;
   FilePipe * pPipe;
   UnchangeableString name;
   char msg[DEF_STR_SZ];
   int i;
   bool found;

   //check adjustable parameters
   for(i = 0; i < m_NumParams; i++)
   {
      name = m_pList[i]->GetName();
      pCur = pList;
      found = false;
      while(pCur != NULL)
      {
         pPipe = pCur->GetPipe();
         if(pPipe->FindAndReplace(name, "0.00") > 0)
         {
            found = true;
            break;
         }         
         pCur = pCur->GetNext();
      }
      if(found == false)
      {
         sprintf(msg, "Parameter |%s| not found in any template file", name);
         LogError(ERR_FILE_IO, msg);
      }
   }

   //this will cause reset of replacement string....
   pCur = pList;
   while(pCur != NULL)
   {
      pPipe = pCur->GetPipe();
      pPipe->StringToFile();
      pCur = pCur->GetNext();
   }
}/* end CheckTemplateFiles() */

/******************************************************************************
CheckMnemonics()

Checks to see if and parameter name is nested within another parameter name 
(e.g. Kback is nested within Kbackground). Since the parameter substitution 
routine uses strstr() such nesting can cause undesirable behavior and should 
be reported as an error to the user.
******************************************************************************/
void ParameterGroupWorker::CheckMnemonics(void)
{
   UnchangeableString name;
   UnchangeableString comp;
   char msg[DEF_STR_SZ];
   int i, j;
   bool found = false;

   //check adjustable parameters
   for(i = 0; i < m_NumParams; i++)
   {
      name = m_pList[i]->GetName();      

      //check adjustable parameters
      for(j = 0; j < m_NumParams; j++)
      {
         comp = m_pList[j]->GetName();
         if((i != j) && (strstr(comp, name) != NULL))
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }

   }/* end check adjustable parameters */

   if(found == true)
   {
      ExitProgram(1);
   }
}/* end CheckMnemonics() */

/******************************************************************************
WriteParams()

Stuffs current parameter values using the provided array values. This function
should usually be followed by a call to Model::Execute() to ensure that model
output is consistent with model parameters.

******************************************************************************/
void ParameterGroupWorker::WriteParams(Ironclad1DArray p) {
    double viol = 0.00;
    for (int j = 0; j < m_NumParams; j++) {
        m_pList[j]->SetValue(p[j]);
    }
}/* end WriteParams() */
