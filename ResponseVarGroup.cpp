/******************************************************************************
File      : ResponseVarGroup.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates the response variable group, the group of response variables 
which the objective function (and possibly constraints) is based upon. Response
variables are to optimization, what observations are in regression/calibration.

Version History
05-10-04    lsm   created
01-11-05    lsm   added support for tied response variables
01-20-05    lsm   added support for weighted sum tied resp. vars.
01-01-07    lsm   Revised ValueExtractor CTOR interface. 
******************************************************************************/
#include <string.h>
#include <stdlib.h>

#include "ResponseVarGroup.h"
#include "TiedRespVar.h"
#include "ResponseVar.h"
#include "ValueExtractor.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
GetRespVarPtr()

Returns a pointer to the ith response variable, or NULL if it is out of bounds.
******************************************************************************/
RespVarABC * ResponseVarGroup::GetRespVarPtr(int i)
{  
   if ( i < m_NumRespVars) { return m_pRespVarList[i];}   
   return NULL;
} /* end GetRespVarPtr() */

/******************************************************************************
GetRespVarPtr()

Returns a pointer to the  response variable matching the name arg or NULL if 
no match found.
******************************************************************************/
RespVarABC * ResponseVarGroup::GetRespVarPtr(IroncladString name)
{  
   int i;

   for(i = 0; i < m_NumRespVars; i++)
   {
      if(strcmp(name, m_pRespVarList[i]->GetName()) == 0)
      {
         return m_pRespVarList[i];
      } /* end if() */
   }/* end for() */

   for(i = 0; i < m_NumTiedRespVars; i++)
   {
      if(strcmp(name, m_pTiedRespVarList[i]->GetName()) == 0)
      {
         return m_pTiedRespVarList[i];
      } /* end if() */
   }/* end for() */
  return NULL;
} /* end GetRespVarPtr() */

/******************************************************************************
GetNumRespVars()

Returns the number of  response variables.
******************************************************************************/
int ResponseVarGroup::GetNumRespVars(void)
{
  return m_NumRespVars;
}/* end GetNumRespVars() */

/******************************************************************************
GetNumTiedRespVars()

Returns the number of tied response variables.
******************************************************************************/
int ResponseVarGroup::GetNumTiedRespVars(void)
{
  return m_NumTiedRespVars;
}/* end GetNumTiedRespVars() */

/******************************************************************************
WriteList()

Writes the details of all the response variables.
******************************************************************************/
void ResponseVarGroup::WriteList(FILE * pFile, int type)
{
   int i;

   for(i = 0; i < m_NumRespVars; i++)
   {     
      m_pRespVarList[i]->Write(pFile, type);     
   }/* end for() */

   for(i = 0; i < m_NumTiedRespVars; i++)
   {     
      m_pTiedRespVarList[i]->Write(pFile, type);     
   }/* end for() */
} /* end WriteList() */

/******************************************************************************
ExtractVals()

Extracts values for each response variable from the corresponding output file
and stores them in the current value field.
******************************************************************************/
void ResponseVarGroup::ExtractVals(void)
{ 
   char errMsg[DEF_STR_SZ];
   bool bOk;
   int i;
   double currentValue;
   UnchangeableString name;
   int line;
   int col;
   char tok;
   UnchangeableString keyword;   

   //read output files into memory
   m_pRespFiles->ReadOutputFiles();
   
   for(i = 0; i < m_NumRespVars; i++)
   {     
      name = m_pRespVarList[i]->GetFileName();
      line = m_pRespVarList[i]->GetLine();
      col = m_pRespVarList[i]->GetColumn();
      keyword = m_pRespVarList[i]->GetKeyword();
      tok = m_pRespVarList[i]->GetToken();
      bOk = m_pRespFiles->ExtractValue(name, keyword, line, col, tok, &currentValue);
	  if(bOk == false)
	  {
		 LogError(ERR_CONTINUE,"Ostrich failed to process the following response variable:");
		 sprintf(errMsg, "Name    : %s", m_pRespVarList[i]->GetName());
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "File    : %s", name);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Line    : %d", line);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Column  : %d", col);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Keyword : %s", keyword);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Token   : '%c'", tok);
		 LogError(ERR_CONTINUE, errMsg);
	     ExitProgram(1);
	  }
      m_pRespVarList[i]->SetCurrentVal(currentValue);
   } /* end for() */
} /* end ExtractVals() */

/******************************************************************************
InitializeVals()

Extracts values for each response variable from the corresponding output file and
stores them in the initial value field.
******************************************************************************/
void ResponseVarGroup::InitializeVals(void)
{ 
   char errMsg[DEF_STR_SZ];
   bool bOk;
   int i;
   double initialValue;
   UnchangeableString name;
   int line;
   int col;
   char tok;
   UnchangeableString keyword;

   //read output files into memory
   m_pRespFiles->ReadOutputFiles();
   
   for(i = 0; i < m_NumRespVars; i++)
   {     
      name = m_pRespVarList[i]->GetFileName();
      line = m_pRespVarList[i]->GetLine();
      col = m_pRespVarList[i]->GetColumn();
      keyword = m_pRespVarList[i]->GetKeyword();
      tok = m_pRespVarList[i]->GetToken();
      bOk = m_pRespFiles->ExtractValue(name, keyword, line, col, tok, &initialValue);
	  if(bOk == false)
	  {
		 LogError(ERR_CONTINUE,"Ostrich failed to process the following response variable:");
		 sprintf(errMsg, "Name    : %s", m_pRespVarList[i]->GetName());
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "File    : %s", name);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Line    : %d", line);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Column  : %d", col);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Keyword : %s", keyword);
		 LogError(ERR_CONTINUE, errMsg);
		 sprintf(errMsg, "Token   : '%c'", tok);
		 LogError(ERR_CONTINUE, errMsg);
	     ExitProgram(1);
	  }
      m_pRespVarList[i]->SetInitialVal(initialValue);
   } /* end for() */
} /* end InitializeVals() */

/******************************************************************************
CTOR

Associates the object with an input file containing the details of each 
response variable.
******************************************************************************/
ResponseVarGroup::ResponseVarGroup(void)
{
   m_NumRespVars = 0;
   m_pRespVarList = NULL;
   m_NumTiedRespVars = 0;
   m_pTiedRespVarList = NULL;
   m_pRespFiles = NULL;   

   InitFromFile(GetInFileName());

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
CTOR

Associates the object with an input file containing the details of each 
response variable. Overrides the default Response Variable section tokens.
******************************************************************************/
ResponseVarGroup::ResponseVarGroup(const char * pToken)
{
   char * pStart;
   char * pEnd;
   pStart = new char[strlen("Begin")+strlen(pToken)+1]; pStart[0] = (char)NULL;
   pEnd   = new char[strlen("End")+strlen(pToken)+1]; pEnd[0] = (char)NULL;
   strcpy(pStart, "Begin"); strcat(pStart, pToken);
   strcpy(pEnd, "End"); strcat(pEnd, pToken);
   
   m_NumRespVars = 0;
   m_pRespVarList = NULL;
   m_NumTiedRespVars = 0;
   m_pTiedRespVarList = NULL;
   m_pRespFiles = NULL;   

   InitFromFile(GetInFileName(), pStart, pEnd);
   delete [] pStart;
   delete [] pEnd;

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()
******************************************************************************/   
void ResponseVarGroup::Destroy(void)
{
   int i;
    
   for(i = 0; i < m_NumRespVars; i++)
   {
      delete m_pRespVarList[i];
   }
   delete [] m_pRespVarList;

   for(i = 0; i < m_NumTiedRespVars; i++) 
   {
      delete m_pTiedRespVarList[i];
   }
   delete [] m_pTiedRespVarList;

   delete m_pRespFiles;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
InitFromFile()

Reads the response variable data from the input file.
******************************************************************************/
void ResponseVarGroup::InitFromFile(IroncladString respFileName)
{
   const char * start_tag = "BeginResponseVars";
   const char * end_tag   = "EndResponseVars";
   int    start_len = (int)strlen(start_tag);
   int    end_len   = (int)strlen(end_tag);

   int i, j; 
   char tmpName[DEF_STR_SZ];
   char tmpKey[DEF_STR_SZ];
   int line;
   int col;
   char tmpFile[DEF_STR_SZ];
   char * lineStr;
   char * pTok, tok;
   char tmp1[DEF_STR_SZ];
   FILE * pFile;
   ParameterABC * pLine, * pCol;
   TiedParamABC * pTiedLine, * pTiedCol;
   bool bAug;

   pFile = fopen(respFileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("ResponseVarGroup::InitFromFile", respFileName);
   }/* end if() */

   //make sure correct tokens are present
   FindToken(pFile, start_tag, respFileName);
   FindToken(pFile, end_tag, respFileName);
   rewind(pFile);

   //count number of response variables
   m_NumRespVars = 0;
   FindToken(pFile, start_tag, respFileName);
   lineStr = GetNxtDataLine(pFile, respFileName);

   while(strncmp(lineStr, end_tag, end_len) != 0)
   {   
      m_NumRespVars++;
      lineStr = GetNxtDataLine(pFile, respFileName);
   }/* end while() */
   rewind(pFile);

   if(m_NumRespVars == 0)
   {
      LogError(ERR_FILE_IO,"No response variables specified");
      return;
   }/* end if() */

   //read in each response variable
   i = 0;
   NEW_PRINT("ResponseVar *", m_NumRespVars);
   m_pRespVarList = new ResponseVar * [m_NumRespVars];

   FindToken(pFile,start_tag, respFileName);
   lineStr = GetNxtDataLine(pFile, respFileName);
   while(strncmp(lineStr, end_tag, end_len) != 0)
   {      
      pTok = lineStr;
      //extract name of response variable (no spaces allowed)
      j = ExtractString(pTok, tmpName);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitFromFile()");
      pTok += j;
      //extract filename (spaces allowed)
      j = ExtractFileName(pTok, tmpFile);
      pTok += j;
      //extract keyword
      j = ExtractString(pTok, tmpKey);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitFromFile()");
      pTok += j;
      //extract line
      j = ExtractString(pTok, tmp1);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitFromFile()");
      line = atoi(tmp1);
      
      pLine = GetParameterByName(tmp1);
      pTiedLine = GetTiedParameterByName(tmp1);      
      pTok += j;
      //extract column
      j = ExtractString(pTok, tmp1);
      col = atoi(tmp1);
      pCol = GetParameterByName(tmp1);
      pTiedCol = GetTiedParameterByName(tmp1);

      if(j == -1) //no token specified, default to whitespace.
	  {
		  tok = ' ';
	  }
	  else
	  {		
         pTok += j;

		 if (strstr(pTok, "' '") != NULL) // pTok == ' ' (whitespace separated by single quotes)
	     {
		    tok = ' ';
		    pTok = strstr(pTok, "' '");
			pTok += 3;
	     }
	     else //extract token
	     {
		     //extract token
		     j = ExtractString(pTok, tmp1);
		     if((tmp1[0] == 0x27) && (tmp1[2] == 0x27)){ tok = tmp1[1];}  //(wrapped in ' chars)
		     else if (tmp1[0] != (char)NULL) tok = tmp1[0]; //not wrapped in quotes
		     else{ tok = ' ';} //default
		     pTok += j;
	     }

         //extract augmented output flag
         bAug;
         j = ExtractString(pTok, tmp1);
         if(strcmp(tmp1, "yes") == 0){ bAug = true; }
         else{ bAug = false;}
         pTok += j;
	  }
   
      NEW_PRINT("ResponseVar", 1);
      m_pRespVarList[i] = new ResponseVar(tmpName,tmpFile,tmpKey,line,col,tok,bAug);
      m_pRespVarList[i]->SetColPtr(pCol);
      m_pRespVarList[i]->SetLinePtr(pLine);
      m_pRespVarList[i]->SetColPtr(pTiedCol);
      m_pRespVarList[i]->SetLinePtr(pTiedLine);

      /*------------------------------------------------------------------
      Create a ValueExtractor class for the given file (if one with that 
      name (i.e. tmpFile) hasn't already been created).
      ------------------------------------------------------------------*/      
      if(m_pRespFiles == NULL)
      { 
         NEW_PRINT("ValueExtractor", 1);
         m_pRespFiles = new ValueExtractor(tmpFile, true, 0.00);
      }
      else{ m_pRespFiles->Insert(tmpFile);}

      i++;
      
      lineStr = GetNxtDataLine(pFile, respFileName);
   }/* end while() */
   
   fclose(pFile);

   InitTiedRespVars(respFileName);
} /* end InitFromFile() */

/******************************************************************************
InitFromFile()

Reads the response variable data from the input file. Section tokens are passed
in as pStart and pEnd
******************************************************************************/
void ResponseVarGroup::InitFromFile(IroncladString respFileName, const char * start_tag, const char * end_tag)
{
   int    start_len = (int)strlen(start_tag);
   int    end_len   = (int)strlen(end_tag);

   int i, j; 
   char tmpName[DEF_STR_SZ];
   char tmpKey[DEF_STR_SZ];
   int line;
   int col;
   char tmpFile[DEF_STR_SZ];
   char * lineStr;
   char * pTok, tok;
   char tmp1[DEF_STR_SZ];
   FILE * pFile;
   ParameterABC * pLine, * pCol;
   TiedParamABC * pTiedLine, * pTiedCol;
   bool bAug;

   pFile = fopen(respFileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("ResponseVarGroup::InitFromFile", respFileName);
   }/* end if() */

   if(!CheckToken(pFile, start_tag, respFileName))
   {
      fclose(pFile);
      return;
   }
   //make sure correct tokens are present
   rewind(pFile);
   FindToken(pFile, start_tag, respFileName);
   FindToken(pFile, end_tag, respFileName);
   rewind(pFile);

   //count number of response variables
   m_NumRespVars = 0;
   FindToken(pFile, start_tag, respFileName);
   lineStr = GetNxtDataLine(pFile, respFileName);

   while(strncmp(lineStr, end_tag, end_len) != 0)
   {   
      m_NumRespVars++;
      lineStr = GetNxtDataLine(pFile, respFileName);
   }/* end while() */
   rewind(pFile);

   if(m_NumRespVars == 0)
   {
      LogError(ERR_FILE_IO,"No response variables specified");
      return;
   }/* end if() */

   //read in each response variable
   i = 0;
   NEW_PRINT("ResponseVar *", m_NumRespVars);
   m_pRespVarList = new ResponseVar * [m_NumRespVars];

   FindToken(pFile,start_tag, respFileName);
   lineStr = GetNxtDataLine(pFile, respFileName);
   while(strncmp(lineStr, end_tag, end_len) != 0)
   {      
      pTok = lineStr;
      //extract name of response variable (no spaces allowed)
      j = ExtractString(pTok, tmpName);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitFromFile()");
      pTok += j;
      //extract filename (spaces allowed)
      j = ExtractFileName(pTok, tmpFile);
      pTok += j;
      //extract keyword
      j = ExtractString(pTok, tmpKey);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitFromFile()");
      pTok += j;
      //extract line
      j = ExtractString(pTok, tmp1);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitFromFile()");
      line = atoi(tmp1);
      
      pLine = GetParameterByName(tmp1);
      pTiedLine = GetTiedParameterByName(tmp1);      
      pTok += j;
      //extract column
      j = ExtractString(pTok, tmp1);
      col = atoi(tmp1);
      pCol = GetParameterByName(tmp1);
      pTiedCol = GetTiedParameterByName(tmp1);

      if(j == -1) //no token specified, default to whitespace.
	   {
		   tok = ' ';
	   }
	   else
	   {		
         pTok += j;

		   if (strstr(pTok, "' '") != NULL) // pTok == ' ' (whitespace separated by single quotes)
	      {
		      tok = ' ';
		      pTok = strstr(pTok, "' '");
			   pTok += 3;
	      }
	      else //extract token
	      {
		      //extract token
		      j = ExtractString(pTok, tmp1);
		      if((tmp1[0] == 0x27) && (tmp1[2] == 0x27)){ tok = tmp1[1];}  //(wrapped in ' chars)
		      else if (tmp1[0] != (char)NULL) tok = tmp1[0]; //not wrapped in quotes
		      else{ tok = ' ';} //default
		      pTok += j;
	      }

         //extract augmented output flag
         bAug;
         j = ExtractString(pTok, tmp1);
         if(strcmp(tmp1, "yes") == 0){ bAug = true; }
         else{ bAug = false;}
         pTok += j;
	   }/* end else(token) */
   
      NEW_PRINT("ResponseVar", 1);
      m_pRespVarList[i] = new ResponseVar(tmpName,tmpFile,tmpKey,line,col,tok,bAug);
      m_pRespVarList[i]->SetColPtr(pCol);
      m_pRespVarList[i]->SetLinePtr(pLine);
      m_pRespVarList[i]->SetColPtr(pTiedCol);
      m_pRespVarList[i]->SetLinePtr(pTiedLine);

      /*------------------------------------------------------------------
      Create a ValueExtractor class for the given file (if one with that 
      name (i.e. tmpFile) hasn't already been created).
      ------------------------------------------------------------------*/      
      if(m_pRespFiles == NULL)
      { 
         NEW_PRINT("ValueExtractor", 1);
         m_pRespFiles = new ValueExtractor(tmpFile, true, 0.00);
      }
      else{ m_pRespFiles->Insert(tmpFile);}

      i++;
      
      lineStr = GetNxtDataLine(pFile, respFileName);
   }/* end while() */
   
   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
InitTiedRespVars()

Reads the tied response variable data from the input file.
******************************************************************************/
void ResponseVarGroup::InitTiedRespVars(IroncladString pFileName)
{
   const char * start_tag = "BeginTiedRespVars";
   const char * end_tag   = "EndTiedRespVars";
   int    start_len = (int)strlen(start_tag);
   int    end_len   = (int)strlen(end_tag);

   int i, j, n, nrv;
   FILE * pFile;
   char * pTok;
   RespVarABC ** pTies;
   char nameStr[DEF_STR_SZ]; //name of tied parameter
   char typeStr[DEF_STR_SZ]; //type of relationship
   char * lineStr;
   char tmpStr[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   bool invNRV;

   pFile = fopen(pFileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("ResponseVarGroup::InitTiedRespVars", pFileName);
   }/* end if() */

   //check for token
   if(CheckToken(pFile, start_tag, pFileName) == false)
   {
      fclose(pFile);
      return;
   }
   FindToken(pFile, end_tag, pFileName);
   rewind(pFile);

   //count number of tied response variables
   m_NumTiedRespVars = 0;
   FindToken(pFile, start_tag, pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);

   while(strncmp(lineStr, end_tag, end_len) != 0)
   {   
      m_NumTiedRespVars++;
      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */
   rewind(pFile);

   if(m_NumTiedRespVars == 0)
   {
      LogError(ERR_FILE_IO,"No tied response variables specified");
      return;
   }/* end if() */

   //read in each response variable
   i = 0;
   NEW_PRINT("RespVarABC *", m_NumTiedRespVars);
   m_pTiedRespVarList = new RespVarABC * [m_NumTiedRespVars];
   MEM_CHECK(m_pTiedRespVarList);
   for(i = 0; i < m_NumTiedRespVars; i++){ m_pTiedRespVarList[i] = NULL;}

   i = 0;   
   FindToken(pFile, start_tag, pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);
   while(strncmp(lineStr, end_tag, end_len) != 0)
   {      
      pTok = lineStr;
      //extract name (no spaces allowed)
      j = ExtractString(pTok, nameStr);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitTiedRespVars()");
      pTok += j;
      //extract number of resp. vars.
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitTiedRespVars()");
      pTok += j;
      invNRV = false;
      nrv = atoi(tmpStr);
      //extract names of the resp. vars.
      NEW_PRINT("RespVarABC *", nrv);
      pTies = new RespVarABC * [nrv];
      MEM_CHECK(pTies);
      for(n = 0; n < nrv; n++)
      {
         j = ExtractString(pTok, tmpStr);
         j = ValidateExtraction(j, n, nrv, "ResponseVarGroup::InitTiedRespVars()");
         pTok += j;
         pTies[n] = GetRespVarPtr(tmpStr);
         if(pTies[n] == NULL)
         {
            sprintf(msg, "InitTiedRespVars(): unknown response variable |%s|", tmpStr);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);
         }         
      }/* end for() */
      //extract type of relationship
      j = ExtractString(pTok, typeStr);
      j = ValidateExtraction(j, 1, 1, "ResponseVarGroup::InitTiedRespVars()");
      pTok += j;      
                  
      //pass the rest of the line to the appropriate CTOR
      if(strcmp(typeStr, "linear") == 0)
      { 
         if(nrv == 1)
         { 
            NEW_PRINT("TiedRespVarLin1", 1);
            m_pTiedRespVarList[i] = new TiedRespVarLin1(nameStr, pTies[0], pTok);
         }
         else if(nrv == 2)
         { 
            NEW_PRINT("TiedRespVarLin2", 1);
            m_pTiedRespVarList[i] = new TiedRespVarLin2(nameStr, pTies[0], pTies[1], pTok);
         }
         else{ invNRV = true;}
      }
      else if(strcmp(typeStr, "wsum") == 0)
      { 
         NEW_PRINT("TiedRespVarWsum", 1);
         m_pTiedRespVarList[i] = new TiedRespVarWsum(nameStr, pTies, nrv, pTok);
      }
      else
      {
         sprintf(msg, "InitTiedRespVars(): unknown relationship type |%s|", typeStr);
         LogError(ERR_FILE_IO, msg);
         ExitProgram(1);
      }

      /*
      Make sure number of parameters is compatible with existing 
      tied parameter relationships
      */
      if(invNRV == true)
      {
         sprintf(msg, "InitTiedRespVars(): invalid # of response variables (%d) for type (%s)", nrv, typeStr);
         LogError(ERR_FILE_IO, msg);
         ExitProgram(1);
      }

      MEM_CHECK(m_pTiedRespVarList[i]);

      delete [] pTies;
      lineStr = GetNxtDataLine(pFile, pFileName);
      i++;
   }/* end while() */

   fclose(pFile);
}/* end InitTiedRespVars() */

/******************************************************************************
Write()

Writes user-specified simulated output to the pFile argument.
******************************************************************************/
void ResponseVarGroup::Write(FILE * pFile, int type)
{
   int i;
   
   for(i = 0; i < m_NumRespVars; i++)
   {
      if(m_pRespVarList[i]->IsAugmented())
      {
         m_pRespVarList[i]->WriteSim(pFile, type);
      }
   }
} /* end write() */
