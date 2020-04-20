/******************************************************************************
File      : SurrogateParameterGroup.cpp
Author    : L. Shawn Matott
Copyright : 2006, L. Shawn Matott

Encapsulates a group of surrogate parameters. These are parameters that are
tied to the parameters of the complex model.

Version History
04-06-06    lsm   added copyright information and initial comments.
******************************************************************************/
#include <string.h>
#include <stdlib.h>

#include "SurrogateParameterGroup.h"
#include "ParameterGroup.h"
#include "TiedParamABC.h"
#include "FilePipe.h"
#include "FilePair.h"

#include "Utility.h"
#include "Exception.h"

/******************************************************************************
GetTiedParamPtr()

Retrieves a pointer to the tied parameter with matching name.
******************************************************************************/
TiedParamABC * SurrogateParameterGroup::GetTiedParamPtr(IroncladString name)
{
   int j;

   //determine indices by examing parameter names
   for(j = 0; j < m_NumTied; j++)
   {
      if(strcmp(m_pTied[j]->GetName(), name) == 0){ return m_pTied[j];}
   }/* end for() */
   return NULL;
} /* end GetTiedParamPtr() */

/******************************************************************************
Destroy()

Frees the memory used by the objects of all the parameters contained in it
******************************************************************************/
void SurrogateParameterGroup::Destroy(void)
{
   int j;

   for(j = 0; j < m_NumTied; j++)
   {
      delete m_pTied[j];
   }
   delete [] m_pTied;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR

Initializes parameter group from user-specified input file.
******************************************************************************/
SurrogateParameterGroup::SurrogateParameterGroup
(
   UnmoveableString pFileName, 
   ParameterGroup * pComplex
)
{
   m_pTied = NULL;
   m_NumTied = 0;
   InitTiedParams(pFileName, pComplex);
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SubIntoFile()

Substitutes the estimated value of the parameter into the model input file.
******************************************************************************/
void SurrogateParameterGroup::SubIntoFile(FilePipe * pPipe)
{ 
   int i;
   char find[DEF_STR_SZ];
   char replace[DEF_STR_SZ];
   TiedParamABC * pTied;

   //Tied parameters
   for(i = 0; i < m_NumTied; i++)
   {
      pTied = m_pTied[i];
      strcpy(find, pTied->GetName());
      pTied->GetValAsStr(replace);
      pPipe->FindAndReplace(find,replace);
   } /* end for() */

   pPipe->StringToFile();
} /* end SubIntoFile() */

/******************************************************************************
InitTiedParams()

Reads tied parameter detail from a file.
******************************************************************************/
void SurrogateParameterGroup::InitTiedParams(IroncladString pFileName, ParameterGroup * pComplex)
{
   int i, j, n, np;
   FILE * pFile;
   char * pTok;
   MetaParameter* pParams;
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
   FindToken(pFile, "BeginTiedParams", pFileName);
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

   //abort if no parameters in the section
   if(m_NumTied == 0)
   {
      LogError(ERR_IN_PARSE, "Surrogate model must have at least one tied parameter");
      fclose(pFile);
      ExitProgram(1);
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
      j = ValidateExtraction(j, 1, 1, "SurrogateParameterGroup::InitTiedParams()");
      pTok += j;
      //extract number of parameters
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, 1, 1, "SurrogateParameterGroup::InitTiedParams()");
      pTok += j;
      invalidNumParams = false;
      np = atoi(tmpStr);
      //extract names of the parameters
      NEW_PRINT("MetaParameter", np);
      pParams = new MetaParameter[np];
      for(n = 0; n < np; n++)
      {
         j = ExtractString(pTok, tmpStr);
         j = ValidateExtraction(j, n, np, "SurrogateParameterGroup::InitTiedParams()");
         pTok += j;
         pParams[n] = pComplex->GetMetaParam(tmpStr);
         if(pParams[n].pParam == NULL)
         {
            sprintf(msg, "InitTiedParams(): unknown parameter |%s|", tmpStr);
            LogError(ERR_FILE_IO, msg);
            ExitProgram(1);
         }         
      }/* end for() */
      //extract type of relationship
      j = ExtractString(pTok, typeStr);
      j = ValidateExtraction(j, 1, 1, "SurrogateParameterGroup::InitTiedParams()");
      pTok += j;      
                  
      //pass the entire parameter line to the appropriate CTOR
      if(strcmp(typeStr, "linear") == 0)
      { 
         if(np == 1)
         { 
            NEW_PRINT("TiedParamLin1", 1);
            m_pTied[i] = new TiedParamLin1(nameStr, &pParams[0], pTok);
         }
         else if(np == 2)
         { 
            NEW_PRINT("TiedParamLin2", 1);
            m_pTied[i] = new TiedParamLin2(nameStr, &pParams[0], &pParams[1], pTok);
         }
         else{ invalidNumParams = true;}
      }
	  else if(strcmp(typeStr, "ratio") == 0)
      { 
         if(np == 2)
         { 
            NEW_PRINT("TiedParamSimpleRatio", 1);
            m_pTied[i] = new TiedParamSimpleRatio(nameStr, &pParams[0], &pParams[1], pTok);
         }
         else if(np == 3)
         {
            NEW_PRINT("TiedParamComplexRatio", 1);
            m_pTied[i] = new TiedParamComplexRatio(nameStr, &pParams[0], &pParams[1], &pParams[2], pTok);
         }
         else{ invalidNumParams = true;}
      }
      else if(strcmp(typeStr, "exp") == 0)
      { 
         if(np == 1)
         { 
            NEW_PRINT("TiedParamExp", 1);
            m_pTied[i] = new TiedParamExp(nameStr, &pParams[0], pTok);
         }
         else{ invalidNumParams = true;}
      }
      else if(strcmp(typeStr, "log") == 0)
      { 
         if(np == 1)
         { 
            NEW_PRINT("TiedParamLog", 1);
            m_pTied[i] = new TiedParamLog(nameStr, &pParams[0], pTok);
         }
         else{ invalidNumParams = true;}
      }
      else if(strcmp(typeStr, "dist") == 0)
      { 
         if(np == 4)
         { 
            NEW_PRINT("TiedDistXY", 1);
            m_pTied[i] = new TiedDistXY(nameStr, &pParams[0], &pParams[1], &pParams[2], &pParams[3], pTok);
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

      MEM_CHECK(m_pTied[i]);

      delete [] pParams;
      lineStr = GetNxtDataLine(pFile, pFileName);   
      i++;
   }/* end while() */

   fclose(pFile);
} /* end InitTiedParams() */

/******************************************************************************
Write()

Writes formatted output to the pFile argument.
******************************************************************************/
void SurrogateParameterGroup::Write(FILE * pFile, int type)
{
   int i;
   
   for(i = 0; i < m_NumTied; i++)
   {      
      m_pTied[i]->Write(pFile, type);
   }
} /* end writeToFile() */

/******************************************************************************
CheckTemplateFiles()

Checks to see if every parameter is included in at least one template file.
Parameters not found in any template file will trigger a warning message but 
will not halt the program.
******************************************************************************/
void SurrogateParameterGroup::CheckTemplateFiles(FilePair * pList)
{
   FilePair * pCur;
   FilePipe * pPipe;
   UnchangeableString name;
   char msg[DEF_STR_SZ];
   int i;
   bool found;

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

   //this will cause reset of replacement string....
   pCur = pList;
   while(pCur != NULL)
   {
      pPipe = pCur->GetPipe();
      pPipe->StringToFile();
      pCur = pCur->GetNext();
   }
}/* end CheckTemplateFiles() */

