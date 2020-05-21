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

#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "TiedParamABC.h"
#include "GeomParamABC.h"
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
int ParameterGroup::GetNumParams(void)
{
  return m_NumParams;
} /* end GetNumParams() */

/******************************************************************************
ReadParams()

Stuffs an array with the current parameter values. Array must have been 
previously allocated.
******************************************************************************/
void ParameterGroup::ReadParams(double * p)
{
   for(int j = 0; j < m_NumParams; j++)
   {
      p[j] = m_pList[j]->GetEstVal();
   }
}/* end ReadParams() */

/******************************************************************************
WriteParams()

Stuffs current parameter values using the provided array values. This function
should usually be followed by a call to Model::Execute() to ensure that model 
output is consistent with model parameters.

Returns the amount of violation of parameter bounds (if p contains values that
are outside the parameter limits).
******************************************************************************/
double ParameterGroup::WriteParams(Ironclad1DArray p)
{
   double viol = 0.00;
   for(int j = 0; j < m_NumParams; j++)
   {
      viol += m_pList[j]->SetEstVal(p[j]); 
   }
   return viol;
}/* end WriteParams() */

/******************************************************************************
GetParamPtr()

Retrieves a pointer to the ith parameter.
******************************************************************************/
ParameterABC * ParameterGroup::GetParamPtr(int i)
{
  return m_pList[i];
} /* end GetParamPtr() */

/******************************************************************************
GetParamPtr()

Retrieves a pointer to the parameter with matching name.
******************************************************************************/
ParameterABC * ParameterGroup::GetParamPtr(IroncladString name)
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
GetMetaParam()

Retrieves a parameter (regular or tied) with matching name, packed as a meta-
parameter, meaning both pointer and type.
******************************************************************************/
MetaParameter ParameterGroup::GetMetaParam(IroncladString name)
{
   MetaParameter mp;
   mp.type = BAD_PARAMETER;

   mp.pParam = (void *)(GetParamPtr(name));
   if(mp.pParam != NULL)
   {
      mp.type = RGLR_PARAMETER;
      return mp;
   }

   mp.pParam = (void *)(GetTiedParamPtr(name));
   if(mp.pParam != NULL)
   {
      mp.type = TIED_PARAMETER;
      return mp;
   }   

   return mp;
} /* end GetMetaParam() */

/******************************************************************************
GetTiedParamPtr()

Retrieves a pointer to the tied parameter with matching name.
******************************************************************************/
TiedParamABC * ParameterGroup::GetTiedParamPtr(IroncladString name)
{
   int j;

   //determine indices by examing parameter names
   for(j = 0; j < m_NumTied; j++)
   {
       if(m_pTied[j] != NULL)
       {
         if(strcmp(m_pTied[j]->GetName(), name) == 0){ return m_pTied[j];}
       }
   }/* end for() */
   return NULL;
} /* end GetTiedParamPtr() */

/******************************************************************************
ConfigureSpecialParams()

Set the objective function threshold criteria for pre-emption of model
******************************************************************************/
void ParameterGroup::ConfigureSpecialParams(double minObj, double * minCon)
{
   int j;

   for(j = 0; j < m_NumSpecial; j++)
   {
      m_pSpecial[j]->SetEstVal(minObj, minCon[j]);
   }/* end for() */
} /* end ConfigureSpecialParams() */

/******************************************************************************
GetSpecialConstraints()

Retrieve the current values of the response variable associated with each 
constraint. These canbe used in conjunction with ConfigureSpecialParams()
to update pre-emption thresholds.
******************************************************************************/
void ParameterGroup::GetSpecialConstraints(double * pSC)
{
   ConstraintABC * pCon;
   double rv = 0.00;

   int j;

   for(j = 0; j < m_NumSpecial; j++)
   {
      pCon = m_pSpecial[j]->GetConstraint();
	   if(pCon != NULL)
	   {
         rv = pCon->GetResponseVar();
         pSC[j] = rv;
      }
      else
      {
         pSC[j] = 0.00;
      }
   }/* end for() */
}/* end GetSpecialConstraint() */

/******************************************************************************
EnableSpecialParams()

Enable model pre-emption.
******************************************************************************/
void ParameterGroup::EnableSpecialParams(void)
{
	int j;

   for(j = 0; j < m_NumSpecial; j++)
   {
      m_pSpecial[j]->Enable();
   }/* end for() */
}/* end EnableSpecialParams() */

/******************************************************************************
InitSpecialParams()

Example syntax:
	BeginSpecialParams
		#template   initial    special       upper or  cons-
		#mnemonic   value     parameter      lower?    traint
		OST_COST    1E60       BestCost	     n/a      n/a      
		OST_MASS    1E60     BestConstraint  upper    MyPen    
	EndSpecialParams
******************************************************************************/
void ParameterGroup::InitSpecialParams(IroncladString pFileName)
{
   int i;   
   FILE * pFile;
   char pName[DEF_STR_SZ];
   double initialValue;
   char pType[DEF_STR_SZ];
   char pLimit[DEF_STR_SZ];
   char pConstraint[DEF_STR_SZ];
   char * line;
      
   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("InitSpecialParams()", pFileName);}

   //check for parameter token
   if(CheckToken(pFile, "BeginSpecialParams", pFileName) == false)
   { 
      fclose(pFile);
      return;
   }
   FindToken(pFile, "EndSpecialParams", pFileName);
   rewind(pFile);
   
   //count number of parameters
   m_NumSpecial = 0;
   FindToken(pFile, "BeginSpecialParams", pFileName);
   line = GetNxtDataLine(pFile, pFileName);
   while(strstr(line, "EndSpecialParams") == NULL)
   {            
      m_NumSpecial++;
      line = GetNxtDataLine(pFile, pFileName);
   }/* end while() */

   NEW_PRINT("SpecialParam *", m_NumSpecial);
   m_pSpecial = new SpecialParam *[m_NumSpecial];
   MEM_CHECK(m_pSpecial);

   rewind(pFile);
   FindToken(pFile, "BeginSpecialParams", pFileName);
   line = GetNxtDataLine(pFile, pFileName);
   i = 0;
   while(strstr(line, "EndSpecialParams") == NULL)
   {            
      sscanf(line, "%s %lf %s %s %s", pName, &initialValue, pType, pLimit, pConstraint);

      NEW_PRINT("SpecialParam", 1);
      m_pSpecial[i] = new SpecialParam(pName, pType, pLimit, pConstraint, initialValue);
      MEM_CHECK(m_pSpecial[i]);

      line = GetNxtDataLine(pFile, pFileName);
	  i++;
   }/* end while() */

   fclose(pFile);
} /* end InitSpecialParams() */

/******************************************************************************
Destroy()

Frees the memory used by the objects of all the parameters contained in it
******************************************************************************/
void ParameterGroup::Destroy(void)
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

   for(j = 0; j < m_NumTied; j++)
   {
      delete m_pTied[j];
   }
   delete [] m_pTied;

   for(j = 0; j < m_NumGeom; j++)
   {
      delete m_pGeom[j];
   }
   delete [] m_pGeom;

   for(j = 0; j < m_NumSpecial; j++)
   {
      delete m_pSpecial[j];
   }
   delete [] m_pSpecial;

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
ParameterGroup::ParameterGroup(void)
{
   m_pList = NULL;
   m_pExcl = NULL;
   m_pTied = NULL;
   m_pGeom = NULL;
   m_pSpecial = NULL;
   m_ParamNameList = NULL;
   m_NumParams = 0;
   m_NumExcl = 0;
   m_NumTied = 0;
   m_NumGeom = 0;
   m_NumSpecial = 0;
   InitFromFile(GetInFileName());
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SubIntoFile()

Substitutes the estimated value of the parameter into the model input file.
******************************************************************************/
void ParameterGroup::SubIntoFile(FilePipe * pPipe)
{ 
   int i, size;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   char * pRep;
   ParameterABC * pParam;
   TiedParamABC * pTied;
   GeomParamABC * pGeom;
   SpecialParam * pSpecial;

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

   //Tied parameters
   for(i = 0; i < m_NumTied; i++)
   {
      pTied = m_pTied[i];
      strcpy(find, pTied->GetName());
      pTied->GetValAsStr(replace);
      pPipe->FindAndReplace(find,replace);
   } /* end for() */

   //geometry parameters
   for(i = 0; i < m_NumGeom; i++)
   {     
      pGeom = m_pGeom[i];
      strcpy(find, pGeom->GetName());
      size = m_pGeom[i]->GetValStrSize();

      NEW_PRINT("char", size);
      pRep = new char[size+10];
      MEM_CHECK(pRep);

      pGeom->GetValAsStr(pRep);
      pPipe->FindAndReplace(find,pRep);
      delete [] pRep;
   } /* end for() */

   //special parameters
   for(i = 0; i < m_NumSpecial; i++)
   {     
      pSpecial = m_pSpecial[i];
      strcpy(find, pSpecial->GetName());
      pSpecial->GetValAsStr(replace);
      pPipe->FindAndReplace(find,replace);
   } /* end for() */

   pPipe->StringToFile();
} /* end SubIntoFile() */

/******************************************************************************
WriteDatabaseParameter()

Loop over database entries until the desired parameter is written.
******************************************************************************/
void ParameterGroup::WriteDatabaseParameter(DatabaseABC * pDbase, char * find, char * replace)
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
void ParameterGroup::SubIntoDbase(DatabaseABC * pDbase)
{ 
   int i, size;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   char * pRep;
   ParameterABC * pParam;
   TiedParamABC * pTied;
   GeomParamABC * pGeom;
   SpecialParam * pSpecial;

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

   //Tied parameters
   for(i = 0; i < m_NumTied; i++)
   {
      pTied = m_pTied[i];
      strcpy(find, pTied->GetName());
      pTied->GetValAsStr(replace);
      WriteDatabaseParameter(pDbase, find, replace);
   } /* end for() */

   //geometry parameters
   for(i = 0; i < m_NumGeom; i++)
   {     
      pGeom = m_pGeom[i];
      strcpy(find, pGeom->GetName());
      size = m_pGeom[i]->GetValStrSize();

      NEW_PRINT("char", size);
      pRep = new char[size+10];
      MEM_CHECK(pRep);

      pGeom->GetValAsStr(pRep);
      WriteDatabaseParameter(pDbase, find, replace);
      delete [] pRep;
   } /* end for() */

   //special parameters
   for(i = 0; i < m_NumSpecial; i++)
   {     
      pSpecial = m_pSpecial[i];
      strcpy(find, pSpecial->GetName());
      pSpecial->GetValAsStr(replace);
      WriteDatabaseParameter(pDbase, find, replace);
   } /* end for() */
} /* end SubIntoDbase() */

/******************************************************************************
WriteSuperMuseArgs()

This routine is similar in purpose to SubIntoFile(), except that file substution
information (i.e. find/replace pairs) are written to a SuperMUSE arguments file.
Final substitution into the model template files will be performed by a SuperMUSE
client-side tasker script/batch file.
******************************************************************************/
void ParameterGroup::WriteSuperMuseArgs(FILE * pFile)
{
   int i, size;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   char * pRep;
   ParameterABC * pParam;
   TiedParamABC * pTied;
   GeomParamABC * pGeom;

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

   //Tied parameters
   for(i = 0; i < m_NumTied; i++)
   {
      pTied = m_pTied[i];
      strcpy(find, pTied->GetName());
      pTied->GetValAsStr(replace);
      fprintf(pFile, "%s %s ", find, replace);
   } /* end for() */

   //geometry parameters
   for(i = 0; i < m_NumGeom; i++)
   {     
      pGeom = m_pGeom[i];
      strcpy(find, pGeom->GetName());
      size = m_pGeom[i]->GetValStrSize();

      NEW_PRINT("char", size);
      pRep = new char[size+10];
      MEM_CHECK(pRep);

      pGeom->GetValAsStr(pRep);
      fprintf(pFile, "%s %s ", find, pRep);
      delete [] pRep;
   } /* end for() */

   fprintf(pFile, " \n");
}/* end WriteSuperMuseArgs() */

/******************************************************************************
InitFromFile()

Reads parameter details from a file.
******************************************************************************/
void ParameterGroup::InitFromFile(IroncladString pFileName)
{
   //size the parameter list
   m_NumParams = CountParams(pFileName);
   if(m_NumParams == 0)
   {
      LogError(ERR_FILE_IO,"No parameters specified");
      ExitProgram(1);
   }/* end if() */   

   //collect the names so they can be protected during extraction
   m_ParamNameList = new char*[m_NumParams];
   for(int i = 0; i < m_NumParams; i++)
   {
      m_ParamNameList[i] = NULL;
   }
   GetParameterNames(pFileName);

   NEW_PRINT("ParameterABC *", m_NumParams);
   m_pList = new ParameterABC *[m_NumParams];
   MEM_CHECK(m_pList);

   NEW_PRINT("ParameterABC *", m_NumParams);
   m_pExcl = new ParameterABC *[m_NumParams];
   MEM_CHECK(m_pList);

   for(int i = 0; i < m_NumParams; i++)
   { 
      m_pList[i] = NULL;
      m_pExcl[i] = NULL;
   }
   
   InitRealParams(pFileName);
   InitIntParams(pFileName);
   InitComboParams(pFileName);

   InitTiedParams(pFileName);   
   InitGeomParams(pFileName);
} /* end InitFromFile() */

/******************************************************************************
CountParams()

Counts the number of prameters specified in the input file.
******************************************************************************/
int ParameterGroup::CountParams(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   int count = 0;

   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("CountParams()", pFileName); }

   //check for presence of real parameters
   if(CheckToken(pFile, "BeginParams", pFileName) == true)
   {
      FindToken(pFile, "EndParams", pFileName);
      rewind(pFile);
      //count number of parameters
      FindToken(pFile, "BeginParams", pFileName);      
      line = GetNxtDataLine(pFile, pFileName);   
      while(strstr(line, "EndParams") == NULL)
      {            
         count++;      
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
      rewind(pFile);
   }
   
   //check for presence of integer parameters
   if(CheckToken(pFile, "BeginIntegerParams", pFileName) == true)
   {
      FindToken(pFile, "EndIntegerParams", pFileName);
      rewind(pFile);
      //count number of parameters
      FindToken(pFile, "BeginIntegerParams", pFileName);      
      line = GetNxtDataLine(pFile, pFileName);   
      while(strstr(line, "EndIntegerParams") == NULL)
      {            
         count++;      
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
      rewind(pFile);
   }

   //check for presence of combinatorial parameters
   if(CheckToken(pFile, "BeginCombinatorialParams", pFileName) == true)
   {
      FindToken(pFile, "EndCombinatorialParams", pFileName);
      rewind(pFile);
      //count number of parameters
      FindToken(pFile, "BeginCombinatorialParams", pFileName);      
      line = GetNxtDataLine(pFile, pFileName);   
      while(strstr(line, "EndCombinatorialParams") == NULL)
      {            
         count++;      
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
      rewind(pFile);
   }

   fclose(pFile);

   return count;   
} /* end CountParams() */

/******************************************************************************
GetParameterNames()

Get the names of each parameter specified in the input file.
******************************************************************************/
void ParameterGroup::GetParameterNames(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmpname[DEF_STR_SZ];
   int count = 0;

   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("GetParameterNames()", pFileName); }

   //check for presence of real parameters
   if(CheckToken(pFile, "BeginParams", pFileName) == true)
   {
      FindToken(pFile, "EndParams", pFileName);
      rewind(pFile);
      //count number of parameters
      FindToken(pFile, "BeginParams", pFileName);      
      line = GetNxtDataLine(pFile, pFileName);   
      while(strstr(line, "EndParams") == NULL)
      {            
         ExtractString(line, tmpname);
         m_ParamNameList[count] = new char[strlen(tmpname)+1];
         strcpy(m_ParamNameList[count], tmpname);
         count++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
      rewind(pFile);
   }
   
   //check for presence of integer parameters
   if(CheckToken(pFile, "BeginIntegerParams", pFileName) == true)
   {
      FindToken(pFile, "EndIntegerParams", pFileName);
      rewind(pFile);
      //count number of parameters
      FindToken(pFile, "BeginIntegerParams", pFileName);      
      line = GetNxtDataLine(pFile, pFileName);   
      while(strstr(line, "EndIntegerParams") == NULL)
      {            
         ExtractString(line, tmpname);
         m_ParamNameList[count] = new char[strlen(tmpname)+1];
         strcpy(m_ParamNameList[count], tmpname);
         count++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
      rewind(pFile);
   }

   //check for presence of combinatorial parameters
   if(CheckToken(pFile, "BeginCombinatorialParams", pFileName) == true)
   {
      FindToken(pFile, "EndCombinatorialParams", pFileName);
      rewind(pFile);
      //count number of parameters
      FindToken(pFile, "BeginCombinatorialParams", pFileName);      
      line = GetNxtDataLine(pFile, pFileName);   
      while(strstr(line, "EndCombinatorialParams") == NULL)
      {            
         ExtractString(line, tmpname);
         m_ParamNameList[count] = new char[strlen(tmpname)+1];
         strcpy(m_ParamNameList[count], tmpname);
         count++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
      rewind(pFile);
   }

   fclose(pFile);
} /* end GetParameterNames() */

/******************************************************************************
GetNextEmptyParamIdx()

Finds the next unassigned parameter in the array.
******************************************************************************/
int ParameterGroup::GetNextEmptyParamIdx(void)
{
   for(int i = 0; i < m_NumParams; i++){ if(m_pList[i] == NULL) { return i;}}

   LogError(ERR_ARR_BNDS, "GetNextEmptyParamIdx() : array is filled!");   
   return 0;
}/* end GetNextEmptyParamIdx() */

/******************************************************************************
InitRealParams()

Reads continously varying parameter details from a file.
******************************************************************************/
void ParameterGroup::InitRealParams(IroncladString pFileName)
{
   int i;   
   FILE * pFile;
   char tmpName[DEF_STR_SZ];
   double initialValue;
   double upperBound;
   double lowerBound;
   char tmpTrans1[DEF_STR_SZ];
   char tmpTrans2[DEF_STR_SZ];
   char tmpTrans3[DEF_STR_SZ];
   char tmpInitVal[DEF_STR_SZ];
   char tmpFixFmt[DEF_STR_SZ];
   char logMessage[DEF_STR_SZ];
   char * line;
   bool bFixFmt;
      
   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("InitRealParams()", pFileName);}

   //check for parameter token
   if(CheckToken(pFile, "BeginParams", pFileName) == false)
   { 
      fclose(pFile);
      return;
   }
   FindToken(pFile, "EndParams", pFileName);
   rewind(pFile);
   
   FindToken(pFile, "BeginParams", pFileName);
   line = GetNxtDataLine(pFile, pFileName);
   m_bExtracted = false;
   while(strstr(line, "EndParams") == NULL)
   {            
      strcpy(tmpTrans1, "none");
      strcpy(tmpTrans2, "none");
      strcpy(tmpTrans3, "none");
      strcpy(tmpFixFmt, "free");

      sscanf(line, "%s %s %lf %lf %s %s %s %s", tmpName, tmpInitVal, &lowerBound,
             &upperBound, tmpTrans1, tmpTrans2, tmpTrans3, tmpFixFmt);

      if(strcmp(tmpFixFmt, "free") == 0) bFixFmt = false;
      else bFixFmt = true;

      //assign initial value (possibly random)
      if(strcmp(tmpInitVal, "random") == 0){ initialValue = (((double)MyRand() / (double)MY_RAND_MAX) * (upperBound - lowerBound)) + lowerBound;}
      else if(strcmp(tmpInitVal, "extract") == 0)
      {
         bool bSuccess = ExtractInitialValue(tmpName, bFixFmt, &initialValue); 
         if(bSuccess == false) //revert to random assignment if can't extract
         {
            LogError(ERR_FILE_IO, "Coundn't extract parameter value. Dafaulting to random assignment.");
            initialValue = (((double)MyRand() / (double)MY_RAND_MAX) * (upperBound - lowerBound)) + lowerBound;
         }/* end if() */
         else
         {
            sprintf(logMessage, "extracted %s = %E", tmpName, initialValue);
            //LogError(ERR_FILE_IO, logMessage);
         }
      }/* end else if() */
      else{ initialValue = atof(tmpInitVal); }

      i = GetNextEmptyParamIdx();

      NEW_PRINT("RealParam", 1);
      m_pList[i] = new RealParam(tmpName, initialValue, lowerBound, upperBound,
                                 tmpTrans1, tmpTrans2, tmpTrans3, tmpFixFmt);
      MEM_CHECK(m_pList[i]);

       //assign random values from within transformed space
      if(strcmp(tmpInitVal, "random") == 0)
      { 
         upperBound = m_pList[i]->GetUprBnd();
         lowerBound = m_pList[i]->GetLwrBnd();
         initialValue = (((double)MyRand() / (double)MY_RAND_MAX) * (upperBound - lowerBound)) + lowerBound;
         m_pList[i]->SetEstVal(initialValue);
      }

      //assign extracted values from within transformed space
      if(strcmp(tmpInitVal, "extract") == 0)
      { 
         m_pList[i]->SetEstVal(initialValue);
         m_bExtracted = true;
      }
      line = GetNxtDataLine(pFile, pFileName);
   }/* end while() */

   fclose(pFile);
} /* end InitRealParams() */

/******************************************************************************
ExtractInitialValue()

Reads the parameter value from a model input file. Returns true if successful.
******************************************************************************/
bool ParameterGroup::ExtractInitialValue(char * name, bool bFixFmt, double * pVal)
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
InitIntParams()

Reads the integer parameter details from the file.
******************************************************************************/
void ParameterGroup::InitIntParams(IroncladString pFileName)
{
   int i;   
   FILE * pFile;
   char tmpName[DEF_STR_SZ];
   int initialValue;
   int upperBound;
   int lowerBound;
   char * line;
   char tmpInitVal[DEF_STR_SZ];
      
   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("InitIntParams()", pFileName); }

   //check for parameter token
   if(CheckToken(pFile, "BeginIntegerParams", pFileName) == false)
   {
      fclose(pFile);
      return;
   }
   FindToken(pFile, "EndIntegerParams", pFileName);
   rewind(pFile);

   FindToken(pFile, "BeginIntegerParams", pFileName);
   line = GetNxtDataLine(pFile, pFileName);
   while(strstr(line, "EndIntegerParams") == NULL)
   {            
      sscanf(line, "%s %s %d %d", tmpName, tmpInitVal, &lowerBound, &upperBound);

     //assign initial value (possibly random)
      if(strcmp(tmpInitVal, "random") == 0){ initialValue = (int)(((double)MyRand() / (double)MY_RAND_MAX) * (upperBound - lowerBound)) + lowerBound;}
      else{ initialValue = atoi(tmpInitVal); }

      i = GetNextEmptyParamIdx();

      NEW_PRINT("IntParam", 1);
      m_pList[i] = new IntParam(tmpName, initialValue, lowerBound, upperBound);
      MEM_CHECK(m_pList[i]);
      
      line = GetNxtDataLine(pFile, pFileName);
   }/* end while() */

   fclose(pFile);
} /* end InitIntParams() */

/******************************************************************************
InitComboParams()

Reads the combinatorial parameter details from a file.
******************************************************************************/
void ParameterGroup::InitComboParams(IroncladString pFileName)
{
   int i, j;   
   FILE * pFile;
   char * pTok;
   char nameStr[DEF_STR_SZ]; //name of parameter
   char typeStr[DEF_STR_SZ]; //parameter type (string, real or integer)
   char * lineStr;
      
   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("InitComboParams()", pFileName);}

   //check for parameter token
   if(CheckToken(pFile, "BeginCombinatorialParams", pFileName) == false)
   {
      fclose(pFile);
      return ;
   }
   FindToken(pFile, "EndCombinatorialParams",   pFileName);
   rewind(pFile);
   
   FindToken(pFile, "BeginCombinatorialParams", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);
   while(strstr(lineStr, "EndCombinatorialParams") == NULL)
   {      
      pTok = lineStr;
      //extract name of parameter (no spaces allowed)
      j = ExtractString(pTok, nameStr);
      j = ValidateExtraction(j, 1, 1, "InitComboParams()");
      pTok += j;
      //extract type      
      j = ExtractString(pTok, typeStr);
      j = ValidateExtraction(j, 1, 1, "InitComboParams()");
      pTok += j;
            
      i = GetNextEmptyParamIdx();      
      //pass the entire parameter line to the appropriate CTOR
      if(strcmp(typeStr, "real") == 0)
      { 
         NEW_PRINT("ComboDblParam", 1);
         m_pList[i] = new ComboDblParam(nameStr, pTok);
      }
      else if(strcmp(typeStr, "integer") == 0)
      { 
         NEW_PRINT("ComboIntParam", 1);
         m_pList[i] = new ComboIntParam(nameStr, pTok);
      }
      else if(strcmp(typeStr, "string") == 0)
      { 
         NEW_PRINT("ComboStrParam", 1);
         m_pList[i] = new ComboStrParam(nameStr, pTok);
      }
      else
      {
         sprintf(lineStr, "InitComboParams(): unknown combinatorial type |%s|", typeStr);
         LogError(ERR_FILE_IO, lineStr);
         ExitProgram(1);
      }
      MEM_CHECK(m_pList[i]);

      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */

   fclose(pFile);
} /* end InitComboParams() */

/******************************************************************************
InitTiedParams()

Reads tied parameter detail from a file.
******************************************************************************/
void ParameterGroup::InitTiedParams(IroncladString pFileName)
{
   int i, j, n, np;
   FILE * pFile;
   char * pTok;
   MetaParameter * pParams = NULL;
   char nameStr[DEF_STR_SZ]; //name of tied parameter
   char typeStr[DEF_STR_SZ]; //type of realtionship (linear, exp, or log)
   char * lineStr;
   char tmpStr[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   bool invalidNumParams;
      
   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("InitTiedParams()", pFileName);}

   //check for parameter token
   if(CheckToken(pFile, "BeginTiedParams", pFileName) == false)
   {
      fclose(pFile);
      return ;
   }
   FindToken(pFile, "EndTiedParams",   pFileName);
   rewind(pFile);

   //count number of parameters
   FindToken(pFile, "BeginTiedParams", pFileName);      
   lineStr = GetNxtDataLine(pFile, pFileName);   
   while(strstr(lineStr, "EndTiedParams") == NULL)
   {            
      m_NumTied++;      
      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */
   rewind(pFile);

   //abort if no parmaters in the section
   if(m_NumTied == 0)
   {
      fclose(pFile);
      return;
   }
   //allocate space and initialize
   NEW_PRINT("TiedParamABC *", m_NumTied);
   m_pTied = new TiedParamABC * [m_NumTied];
   MEM_CHECK(m_pTied);
   for(i = 0; i < m_NumTied; i++){ m_pTied[i] = NULL;}

   i = 0;   
   FindToken(pFile, "BeginTiedParams", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);
   while(strstr(lineStr, "EndTiedParams") == NULL)
   {      
      pTok = lineStr;
      //extract name of parameter (no spaces allowed)
      j = ExtractString(pTok, nameStr);
      j = ValidateExtraction(j, 1, 1, "InitTiedParams()");
      pTok += j;
      //extract number of parameters
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, 1, 1, "InitTiedParams()");
      pTok += j;
      invalidNumParams = false;
      np = atoi(tmpStr);
	  if(np != 0)
	  {
		//extract names of the parameters
		NEW_PRINT("MetaParameter", np);
		pParams = new MetaParameter[np];
		for(n = 0; n < np; n++)
		{
		   j = ExtractString(pTok, tmpStr);
		   j = ValidateExtraction(j, n, np, "InitTiedParams()");
           pTok += j;
           pParams[n] = GetMetaParam(tmpStr);
           if(pParams[n].pParam == NULL)
           {
              sprintf(msg, "InitTiedParams(): unknown parameter |%s|", tmpStr);
              LogError(ERR_FILE_IO, msg);
              ExitProgram(1);
           }
         }/* end for() */
         //extract type of relationship
         j = ExtractString(pTok, typeStr);
         if(strcmp(typeStr, "dist") != 0)
           j = ValidateExtraction(j, 1, 1, "InitTiedParams()");
         if(j < 0) pTok = NULL;
         else pTok += j;      
                  
         //pass the entire parameter line to the appropriate CTOR
         if(strcmp(typeStr, "linear") == 0)
         { 
            if(np == 1)
            { 
               NEW_PRINT("TiedParamLin1", 1);
               m_pTied[i] = new TiedParamLin1(nameStr, &(pParams[0]), pTok);
            }
            else if(np == 2)
            { 
               NEW_PRINT("TiedParamLin2", 1);
               m_pTied[i] = new TiedParamLin2(nameStr, &(pParams[0]), &(pParams[1]), pTok);
            }
            else{ invalidNumParams = true;}
         }
         else if(strcmp(typeStr, "wsum") == 0)
         { 
            NEW_PRINT("TiedParamWsum", 1);
            m_pTied[i] = new TiedParamWsum(nameStr, pParams, np, pTok);
         }
         else if(strcmp(typeStr, "ratio") == 0)
         { 
            if(np == 2)
            { 
               NEW_PRINT("TiedParamSimpleRatio", 1);
               m_pTied[i] = new TiedParamSimpleRatio(nameStr, &(pParams[0]), &(pParams[1]), pTok);
            }
            else if(np == 3)
            {
               NEW_PRINT("TiedParamComplexRatio", 1);
               m_pTied[i] = new TiedParamComplexRatio(nameStr, &(pParams[0]), &(pParams[1]), &(pParams[2]), pTok);
            }
            else{ invalidNumParams = true;}
         }
         else if(strcmp(typeStr, "exp") == 0)
         { 
            if(np == 1)
            { 
               NEW_PRINT("TiedParamExp", 1);
               m_pTied[i] = new TiedParamExp(nameStr, &(pParams[0]), pTok);
            }
            else{ invalidNumParams = true;}
         }
         else if(strcmp(typeStr, "log") == 0)
         { 
            if(np == 1)
            { 
               NEW_PRINT("TiedParamLog", 1);
               m_pTied[i] = new TiedParamLog(nameStr, &(pParams[0]), pTok);
            }
            else{ invalidNumParams = true;}
         }
         else if(strcmp(typeStr, "dist") == 0)
         { 
            if(np == 4)
            { 
               NEW_PRINT("TiedDistXY", 1);
               m_pTied[i] = new TiedDistXY(nameStr, &(pParams[0]), &(pParams[1]), &(pParams[2]), &(pParams[3]), pTok);
            }
            else{ invalidNumParams = true;}
         }
         else
         {
            sprintf(msg, "InitTiedParams(): unknown relationship type |%s|", typeStr);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);
         }

         /*
         Make sure number of parameters is compatible with existing 
         tied parameter relationships
         */
         if(invalidNumParams == true)
         {
            sprintf(msg, "InitTiedParams(): invalid # of params (%d) for type (%s)", np, typeStr);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);
         }
	  }
	  else // assign constant value to a parameter
	  {
         NEW_PRINT("TiedParamConstant", 1);
         m_pTied[i] = new TiedParamConstant(nameStr, pTok);
	  }

      MEM_CHECK(m_pTied[i]);

      delete [] pParams;
      lineStr = GetNxtDataLine(pFile, pFileName);
      i++;
   }/* end while() */

   fclose(pFile);
} /* end InitTiedParams() */

/******************************************************************************
InitGeomParams()

Reads geometry parameters from a file. Ostrich will preserve the topology of
these parameters, such that:
1) vertices will be inserted wherever two geometry parameters intersect
2) polygon vertices may be re-ordered to ensure a valid polygon is used
******************************************************************************/
void ParameterGroup::InitGeomParams(IroncladString pFileName)
{
   bool isCirc;
   AugVertexList * pNewVert;
   AugCircle * pNewCirc;
   int i, j, num_starts, num_ends;
   FILE * pFile;
   char * pTok;
   char nameStr[DEF_STR_SZ]; //name of geometry parameter
   char typeStr[DEF_STR_SZ]; //type of geometry (poly2, poly3 or line3)
   char * lineStr; 
   char msg[DEF_STR_SZ];   
   char tmpx[DEF_STR_SZ], tmpy[DEF_STR_SZ], tmpz[DEF_STR_SZ], tmpr[DEF_STR_SZ];
      
   pFile = fopen(pFileName, "r");

   //check that file could be opened
   if(pFile == NULL){ FileOpenFailure("InitGeomParams()", pFileName);}

   //check for parameter token
   if(CheckToken(pFile, "BeginGeomParams", pFileName) == false)
   {
      fclose(pFile);
      return ;
   }
   FindToken(pFile, "EndGeomParams", pFileName);
   rewind(pFile);

   //count number of parameters
   FindToken(pFile, "BeginGeomParams", pFileName);      
   lineStr = GetNxtDataLine(pFile, pFileName);   
   num_starts = num_ends = 0;
   while(strncmp(lineStr, "EndGeomParams", 13) != 0)
   {
      if(strncmp(lineStr, "BeginShape", 10) == 0){ num_starts++;}
      if(strncmp(lineStr, "EndShape",    8) == 0){ num_ends++;}
      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */
   rewind(pFile);

   //check for mismatch between number of start and end tags
   if(num_starts != num_ends)
   {      
      LogError(ERR_FILE_IO, "Mismatch between number of BeginShape and EndShape tags");
      ExitProgram(1);         
   }
   else { m_NumGeom = num_starts;}
   //abort if no parmaters in the section
   if(m_NumGeom == 0)
   {
      fclose(pFile);
      return;
   }
   //allocate space and initialize
   NEW_PRINT("GeomParamABC *", m_NumGeom);
   m_pGeom = new GeomParamABC * [m_NumGeom];
   MEM_CHECK(m_pGeom);
   for(i = 0; i < m_NumGeom; i++){ m_pGeom[i] = NULL;}

   i = 0;   
   FindToken(pFile, "BeginGeomParams", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);   
   while(strncmp(lineStr, "EndGeomParams", 13) != 0)
   {
      if(strcmp(lineStr, "BeginShape") == 0)
      {      
         lineStr = GetNxtDataLine(pFile, pFileName);   

         pTok = lineStr;
         //extract parameter name (no spaces allowed)
         j = ExtractString(pTok, nameStr);
         j = ValidateExtraction(j, 1, 1, "InitGeomParams()");
         pTok += j;
         //extract geometry type
         j = ExtractString(pTok, typeStr);         
         pTok += j;

         //initialize class
         isCirc = false;
         if(strcmp(typeStr, "poly2") == 0)
         { 
            NEW_PRINT("GeomParamPoly2", 1);
            m_pGeom[i] = new GeomParamPoly2(nameStr);
            MEM_CHECK(m_pGeom[i]);
         }
         else if(strcmp(typeStr, "poly3") == 0)
         { 
            NEW_PRINT("GeomParamPoly3", 1);
            m_pGeom[i] = new GeomParamPoly3(nameStr);
            MEM_CHECK(m_pGeom[i]);
         }
         else if(strcmp(typeStr, "line2") == 0)
         { 
            LogError(ERR_FILE_IO, "line2 geometry type not supported");
            ExitProgram(1);
         }
         else if(strcmp(typeStr, "line3") == 0)
         { 
            NEW_PRINT("GeomParamLine3", 1);
            m_pGeom[i] = new GeomParamLine3(nameStr);
            MEM_CHECK(m_pGeom[i]);
         }
         else if(strcmp(typeStr, "circ4") == 0)
         { 
            //set flag, CTOR requires further parsing of syntax
            isCirc = true;
         }
         else
         {
            sprintf(msg, "unknown geomtry type |%s|", typeStr);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);
         }         

         //read in vertices
         lineStr = GetNxtDataLine(pFile, pFileName);   
         while(strncmp(lineStr, "EndShape", 8) != 0)
         {      
            pTok = lineStr;
            //extract values as strings (could be numbers or parameter names)
            tmpx[0] = tmpy[0] = tmpz[0] = tmpr[0] = (char)NULL;
            j = ExtractString(pTok, tmpx);
            j = ValidateExtraction(j, 1, 1, "InitGeomParams()");
            pTok += j;
            j = ExtractString(pTok, tmpy);
            j = ValidateExtraction(j, 1, 1, "InitGeomParams()");
            pTok += j;
            j = ExtractString(pTok, tmpz);
            pTok += j;
            //if circle, also get radius
            if(isCirc == true)
            {
               j = ExtractString(pTok, tmpr);
               pTok += j;
               pNewCirc = InitAugCircle(tmpx, tmpy, tmpz, tmpr);
               NEW_PRINT("GeomParamCirc4", 1);               
               m_pGeom[i] = new GeomParamCirc4(nameStr, pNewCirc);
               MEM_CHECK(m_pGeom[i]);
            }
            else
            {
               pNewVert = InitAugVertex(tmpx, tmpy, tmpz);
               m_pGeom[i]->InsertVertex(pNewVert);
            }

            lineStr = GetNxtDataLine(pFile, pFileName);   
         }/* end while() */
         i++;
      }/* end if() */
      else
      {
         lineStr = GetNxtDataLine(pFile, pFileName);   
      }
   } /* end while() */
   fclose(pFile);
} /* end InitGeomParams() */

/******************************************************************************
Write()

Writes formatted output to the pFile argument.
******************************************************************************/
void ParameterGroup::Write(FILE * pFile, int type)
{
   int i;
   
   for(i = 0; i < m_NumParams; i++)
   {      
      m_pList[i]->Write(pFile, type);
   }

   if((type == WRITE_OPT) || (type == WRITE_DBG))
   {
      for(i = 0; i < m_NumTied; i++)
      {      
         m_pTied[i]->Write(pFile, type);
      }
   }
   if(type == WRITE_DBG)
   {
      for(i = 0; i < m_NumGeom; i++)
      {      
         m_pGeom[i]->Write(pFile, type);
      }
   }
} /* end writeToFile() */

/******************************************************************************
CheckTemplateFiles()

Checks to see if every parameter is included in at least one template file.
Parameters not found in any template file will trigger a warning message but 
will not halt the program.
******************************************************************************/
void ParameterGroup::CheckTemplateFiles(FilePair * pList)
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
   //check tied parameters
   for(i = 0; i < m_NumTied; i++)
   {
      name = m_pTied[i]->GetName();
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

   //check geometry parameters
   for(i = 0; i < m_NumGeom; i++)
   {
      name = m_pGeom[i]->GetName();
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
void ParameterGroup::CheckMnemonics(void)
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

      //check tied parameters
      for(j = 0; j < m_NumTied; j++)
      {
         comp = m_pTied[j]->GetName();
         if(strstr(comp, name) != NULL)
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }

      //check geom. parameters
      for(j = 0; j < m_NumGeom; j++)
      {
         comp = m_pGeom[j]->GetName();
         if(strstr(comp, name) != NULL)
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }
   }/* end check adjustable parameters */

   //check tied parameters
   for(i = 0; i < m_NumTied; i++)
   {
      name = m_pTied[i]->GetName();


      //check adjustable parameters
      for(j = 0; j < m_NumParams; j++)
      {
         comp = m_pList[j]->GetName();
         if(strstr(comp, name) != NULL)
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }

      //check tied parameters
      for(j = 0; j < m_NumTied; j++)
      {
         comp = m_pTied[j]->GetName();
         if((i != j) && (strstr(comp, name) != NULL))
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }

      //check geom. parameters
      for(j = 0; j < m_NumGeom; j++)
      {
         comp = m_pGeom[j]->GetName();
         if(strstr(comp, name) != NULL)
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }
   }/* end check tied parameters */

   //check geometry parameters
   for(i = 0; i < m_NumGeom; i++)
   {
      name = m_pGeom[i]->GetName();

      //check adjustable parameters
      for(j = 0; j < m_NumParams; j++)
      {
         comp = m_pList[j]->GetName();
         if(strstr(comp, name) != NULL)
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }

      //check tied parameters
      for(j = 0; j < m_NumTied; j++)
      {
         comp = m_pTied[j]->GetName();
         if(strstr(comp, name) != NULL)
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }

      //check geom. parameters
      for(j = 0; j < m_NumGeom; j++)
      {
         comp = m_pGeom[j]->GetName();
         if((i != j) && (strstr(comp, name) != NULL))
         {
            sprintf(msg, "|%s| is a substring of |%s|", name, comp);
            LogError(ERR_PRM_NEST, msg);
            found = true;
         }
      }
   }/* end check gemo. parameters */

   if(found == true)
   {
      ExitProgram(1);
   }
}/* end CheckMnemonics() */

/******************************************************************************
FixGeometry()

Corrects the topology of the geometry parameters. Returns true if all problems
are successfuly corrected, false otherwise.
******************************************************************************/
bool ParameterGroup::FixGeometry(void)
{
   char msg[DEF_STR_SZ];
   bool bi, b;
   int i, j;

   //if no geometry parameters, nothing to be done
   if(m_NumGeom == 0){ return true;}

   //convert from augmented to normal geometry
   for(i = 0; i < m_NumGeom; i++){ m_pGeom[i]->Convert();}

   //reorder vertices, if necessary
   b = true;
   for(i = 0; i < m_NumGeom; i++)
   {
      bi = m_pGeom[i]->Reorder();
      if(bi == false)
      { 
         sprintf(msg, "geometry reorder failed |%s|", m_pGeom[i]->GetName());
         LogError(ERR_MISMATCH, msg);
         b = false;
      }
   }/* end for() */

   //insert vertices, if necessary
   for(i = 0; i < m_NumGeom; i++)
   {
      for(j = (i+1); j < m_NumGeom; j++)
      {
         bi = m_pGeom[i]->FixVertices(m_pGeom[j]);
         if(bi == false)
         { 
            sprintf(msg, "fix-vertex geometry failed |%s| and |%s|", 
                    m_pGeom[i]->GetName(), m_pGeom[j]->GetName());
            LogError(ERR_MISMATCH, msg);
            b = false;
         }
      }
   }
   return b; 
}/* end FixGeometry() */

/******************************************************************************
InitAugVertex()
******************************************************************************/
AugVertexList * ParameterGroup::InitAugVertex
(IroncladString xstr, IroncladString ystr, IroncladString zstr)
{
   AugVertexList * pList;

   NEW_PRINT("AugVertexList", 1);
   pList = new AugVertexList;
   MEM_CHECK(pList);

   pList->pNxt = NULL;

   pList->x = atof(xstr);
   pList->y = atof(ystr);
   pList->z = atof(zstr);

   pList->px = GetParamPtr(xstr);
   pList->py = GetParamPtr(ystr);
   pList->pz = GetParamPtr(zstr);

   pList->tx = GetTiedParamPtr(xstr);
   pList->ty = GetTiedParamPtr(ystr);
   pList->tz = GetTiedParamPtr(zstr);

   return pList;
}/* end InitAugVertex() */

/******************************************************************************
InitAugCircle()
******************************************************************************/
AugCircle * ParameterGroup::InitAugCircle
(IroncladString xstr, IroncladString ystr, 
 IroncladString zstr, IroncladString rstr)
{
   AugCircle * pCirc;

   NEW_PRINT("AugCircle", 1);
   pCirc = new AugCircle;
   MEM_CHECK(pCirc);

   pCirc->x = atof(xstr);
   pCirc->y = atof(ystr);
   pCirc->z = atof(zstr);
   pCirc->r = atof(rstr);

   pCirc->px = GetParamPtr(xstr);
   pCirc->py = GetParamPtr(ystr);
   pCirc->pz = GetParamPtr(zstr);
   pCirc->pr = GetParamPtr(rstr);

   pCirc->tx = GetTiedParamPtr(xstr);
   pCirc->ty = GetTiedParamPtr(ystr);
   pCirc->tz = GetTiedParamPtr(zstr);
   pCirc->tr = GetTiedParamPtr(rstr);

   return pCirc;
}/* end InitAugCircle() */

/******************************************************************************
CheckBounds()

For each parameter, check that the upper bound is greater than the lower bound.
******************************************************************************/
void ParameterGroup::CheckBounds(void)
{
   int i;
   char msg[DEF_STR_SZ];
   for(i = 0; i < m_NumParams; i++)
   {     
      if(m_pList[i]->GetUprBnd() < m_pList[i]->GetLwrBnd())
      {
         sprintf(msg, "Parameter (%s) has incorrect bounds (upper bound less than lower bound)\n", m_pList[i]->GetName());
         LogError(ERR_FILE_IO, msg);
         ExitProgram(1);
      }
   }
}/* end CheckBounds() */

/******************************************************************************
ExcludeParam()

Remove the given parameter from the calibration by moving it to the exclusion
list.
******************************************************************************/
void ParameterGroup::ExcludeParam(UnchangeableString prm)
{
   int i, j;
   double val;

   for(i = 0; i < m_NumParams; i++)
   {
      if(strcmp(prm, m_pList[i]->GetName()) == 0) break;
   }

   if(i == m_NumParams) return; //no match

   //fix value at midpoint of range
   val = 0.5 * (m_pList[i]->GetUprBnd() + m_pList[i]->GetLwrBnd());
   m_pList[i]->SetEstVal(val);

   //move parameter to excluded list
   m_pExcl[m_NumExcl] = m_pList[i];
   m_NumExcl++;
   
   //remove parameter from active list
   for(j = i; j < (m_NumParams-1); j++)
   {
      m_pList[j] = m_pList[j+1];
   }/* end for() */
   m_pList[j] = NULL;
   m_NumParams--;
}/* end ExcludeParam() */
