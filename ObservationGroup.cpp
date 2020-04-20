/******************************************************************************
File      : ObservationGroup.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates the observation group, the group of observations which the 
objective function is based upon.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added ObsToken
08-17-04    lsm   added reporting of memory allocations, ValueExtractors are
                  now part of ObservationGroup
10-04-04    lsm   moved observation tokens to individual observations
01-01-07    lsm   Added copy CTOR and Read/Write routines to support Surrogate-
                  model approach. Added ExcludeObs() subroutine to support the
                  "hold" observations functionality.
******************************************************************************/
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ObservationGroup.h"
#include "Observation.h"
#include "ValueExtractor.h"

#include "Utility.h"
#include "Exception.h"

/******************************************************************************
ReadObservations()

Stuffs an array with the current simulated observation values. Array must have 
been previously allocated.
******************************************************************************/
void ObservationGroup::ReadObservations(double * obs)
{
   for(int j = 0; j < m_NumObs; j++)
   {
      obs[j] = m_pObsList[j]->GetComputedVal(false, false);
   }
}/* end ReadObservations() */

/******************************************************************************
WriteObservations()

Stuffs current simulated observation values using the provided array values. 
******************************************************************************/
void ObservationGroup::WriteObservations(Ironclad1DArray obs)
{
   for(int j = 0; j < m_NumObs; j++)
   {
      m_pObsList[j]->SetComputedVal(obs[j]); 
   }   
}/* end WriteObservations() */

/******************************************************************************
GetObsPtr()

Returns a pointer to the ith observation, or NULL if i is out of bounds.
******************************************************************************/
Observation * ObservationGroup::GetObsPtr(int i)
{  
   if ( i < m_NumObs) { return m_pObsList[i];}   
   return NULL;
} /* end GetObsPtr() */

/******************************************************************************
GetObsPtr()

Returns a pointer to the observation matching the name arg or NULL if no match
found.
******************************************************************************/
Observation * ObservationGroup::GetObsPtr(IroncladString name)
{  
   int i;

   for(i = 0; i < m_NumObs; i++)
   {
      if(strcmp(name, m_pObsList[i]->GetName()) == 0)
      {
         return m_pObsList[i];
      } /* end if() */
   }/* end for() */
  return NULL;
} /* end GetObsPtr() */

/******************************************************************************
GetNumObs()

Returns the number of observation points.
******************************************************************************/
int ObservationGroup::GetNumObs(void)
{
  return m_NumObs;
}/* end GetNumObs() */

/******************************************************************************
GetNumGroups()

Returns the number of observation groups.
******************************************************************************/
int ObservationGroup::GetNumGroups(void)
{
   if(m_NumGroups > 0) return m_NumGroups;

   UnchangeableString g1;
   UnchangeableString g2;
   bool bFound;
   m_NumGroups = 0;
   for(int i = 0; i < m_NumObs; i++)
   {
      g1 = m_pObsList[i]->GetGroup();
      bFound = false;
      for(int j = (i-1); j >= 0; j--)
      {         
         g2 = m_pObsList[j]->GetGroup();
         if(strcmp(g1,g2) == 0)
         {
            bFound = true;
            break;
         }
      }/* end for() */
      if(bFound == false)
      {
         m_NumGroups++;
      }
   }/* end for() */

  return m_NumGroups;
}/* end GetNumGroups() */

/******************************************************************************
GetGroup()

Get the name of the ith group.
******************************************************************************/
UnchangeableString ObservationGroup::GetGroup(int whichGroup)
{
   if(m_NumGroups == 1) return m_pObsList[0]->GetGroup();
   if(whichGroup == 0)  return m_pObsList[0]->GetGroup();

   UnchangeableString g1;
   UnchangeableString g2;
   bool bFound;
   int g = 0;
   for(int i = 0; i < m_NumObs; i++)
   {
      g1 = m_pObsList[i]->GetGroup();
      bFound = false;
      for(int j = (i-1); j >= 0; j--)
      {       
         g2 = m_pObsList[j]->GetGroup();
         if(strcmp(g1,g2) == 0)
         {
            bFound = true;
            break;
         }
      }/* end for() */
      if(bFound == false)
      {
         if(whichGroup == g)
         {
            return g1;
         }
         g++;
      }
   }/* end for() */
  return NULL;
}/* end GetGroup() */

/******************************************************************************
WriteList()

Writes the details of all the observation points.
******************************************************************************/
void ObservationGroup::WriteList(FILE * pFile, int type)
{
   int i;

   for(i = 0; i < m_NumObs; i++)
   {     
      m_pObsList[i]->Write(pFile, type);     
   }/* end for() */
} /* end PrintList() */

/******************************************************************************
ExtractVals()

Extracts values for each observation from the corresponding output file.
******************************************************************************/
void ObservationGroup::ExtractVals(void)
{ 
   char errMsg[DEF_STR_SZ];
   bool bOk;
   int i;
   double computedValue;
   UnchangeableString name;
   int line;
   int col;
   char tok;
   UnchangeableString keyword;

   //read output files into memory
   m_pObsFiles->ReadOutputFiles();
   
   for(i = 0; i < m_NumObs; i++)
   {     
      name = m_pObsList[i]->GetFileName();
      line = m_pObsList[i]->GetLine();
      col = m_pObsList[i]->GetColumn();
      keyword = m_pObsList[i]->GetKeyword();
      tok = m_pObsList[i]->GetToken();
      bOk = m_pObsFiles->ExtractValue(name, keyword, line, col, tok, &computedValue);
	  if(bOk == false)
	  {
		 LogError(ERR_CONTINUE,"Ostrich failed to process the following observation:");
		 sprintf(errMsg, "Name    : %s", m_pObsList[i]->GetName());
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
      m_pObsList[i]->SetComputedVal(computedValue);
   } /* end for() */
} /* end ExtractVals() */

/******************************************************************************
CTOR

Associates the object with an input file containing the details of each 
observation.
******************************************************************************/
ObservationGroup::ObservationGroup(void)
{
   m_NumObs = 0;
   m_NumGroups = 0;
   m_pObsList = NULL;
   m_pObsFiles = NULL;   

   InitFromFile(GetInFileName());

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Copy CTOR

First copies all relevant settings from the given ObservationGroup, then reads
in filenames and parsing information from the Surrogate model input file.
******************************************************************************/
ObservationGroup::ObservationGroup
(
   ObservationGroup * pCopy,
   UnmoveableString pFileName
)
{
   int i, j; 
   char tmpName[DEF_STR_SZ];
   char tmpKey[DEF_STR_SZ];
   int line;
   int col;
   char tmpFile[DEF_STR_SZ];
   char * lineStr;
   char * pTok, tok;
   char tmp1[DEF_STR_SZ];
   char group[DEF_STR_SZ];
   FILE * pFile;
   bool bQuitOnError;
   double errorVal;

   m_pObsFiles = NULL;

   /*------------------------------------------------------
   Copy information about names, values and weights from
   the complex model. Assign default values to parsing
   information:
      File Name: OST_OBS_FILE
      Keyword  : <obs_name>
      Line     : 0
      Column   : 2
      Token    : whitespace (' ')
   ------------------------------------------------------*/
   m_NumObs = pCopy->GetNumObs();
   m_NumGroups = pCopy->GetNumGroups();

   NEW_PRINT("Observation *", m_NumObs);
   m_pObsList = new Observation * [m_NumObs];
   MEM_CHECK(m_pObsList);

   for(i = 0; i < m_NumObs; i++)
   {
      NEW_PRINT("Observation", 1);
      m_pObsList[i] = new Observation(pCopy->GetObsPtr(i));
      MEM_CHECK(m_pObsList[i]);
   }

   /*-------------------------------------------------------
   Read in parsing information from input file
   -------------------------------------------------------*/
   pFile = fopen(pFileName, "r");

   if(pFile == NULL)
   {
      FileOpenFailure("ObservationGroup::Copy CTOR", pFileName);
   }/* end if() */

   //make sure correct tokens are present
   FindToken(pFile, "BeginObservations", pFileName);
   FindToken(pFile, "EndObservations", pFileName);
   rewind(pFile);

   //check for optional error-handling instructions
   bQuitOnError = true;
   errorVal = 0.00;
   if(CheckToken(pFile, "OnObsError", pFileName) == true)
   {
      lineStr = GetCurDataLine();
      MyStrLwr(lineStr);
      if(strstr(lineStr, "quit") == NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &errorVal);
         bQuitOnError = false;
      }
   }
   rewind(pFile);

   //read in each observation
   FindToken(pFile, "BeginObservations", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);
   while(strstr(lineStr, "EndObservations") == NULL)
   {      
      pTok = lineStr;
      //extract name of observation (no spaces allowed)
      j = ExtractString(pTok, tmpName);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      pTok += j;
      //extract filename (spaces allowed)
      j = ExtractFileName(pTok, tmpFile);
      pTok += j;
      //extract keyword
      j = ExtractString(pTok, tmpKey);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      pTok += j;
      //extract line
      j = ExtractString(pTok, tmp1);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      line = atoi(tmp1);
      pTok += j;
      //extract column
      j = ExtractString(pTok, tmp1);
      col = atoi(tmp1);
      pTok += j;
      //extract token (wrapped in ' chars)
      j = ExtractString(pTok, tmp1);
      if((tmp1[0] == 0x27) && (tmp1[2] == 0x27)){ tok = tmp1[1];}
      else{ tok = ' ';}
      pTok += j;
      //extract augmented output flag
      bool bAug;
      j = ExtractString(pTok, tmp1);
      if(strcmp(tmp1, "yes") == 0) bAug = true;
      else{ bAug = false;}
      pTok += j;
      //extract observation group
      strcpy(group, "none");
      j = ExtractString(pTok, group);
      pTok += j;

      for(i = 0; i < m_NumObs; i++)
      {
         if(strcmp(tmpName, m_pObsList[i]->GetName()) == 0)
         {
            m_pObsList[i]->Reconfigure(tmpFile, tmpKey, line, col, tok, bAug, group);
            break;
         }
      }
      if(i == m_NumObs)
      {
         sprintf(lineStr, "Unknown observation |%s|, no match in complex model",
                 tmpName);
         LogError(ERR_IN_PARSE, lineStr);
         fclose(pFile);
         ExitProgram(1);
      }

      /*------------------------------------------------------------------
      Create a ValueExtractor class for the given file (if one with that 
      name (i.e. tmpFile) hasn't already been created).
      ------------------------------------------------------------------*/      
      if(m_pObsFiles == NULL)
      { 
         NEW_PRINT("ValueExtractor", 1);
         m_pObsFiles = new ValueExtractor(tmpFile, bQuitOnError, errorVal);
      }
      else{ m_pObsFiles->Insert(tmpFile);}

      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */ 
   fclose(pFile);

   /*---------------------------------------------------------------------
   Check to see if there will be any interpolated observations. If so, an 
   appropriate ValueExtractor needs to be inserted.
   ----------------------------------------------------------------------*/
   for(i = 0; i < m_NumObs; i++)
   {
      if(strcmp(m_pObsList[i]->GetFileName(), OST_OBS_FILE) == 0)
      {
         m_pObsFiles->Insert(OST_OBS_FILE);
         break;
      }
   }/* end for() */

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()
******************************************************************************/   
void ObservationGroup::Destroy(void)
{
   int i;
     
   for(i = 0; i < m_NumObs; i++) {delete m_pObsList[i];}
   delete [] m_pObsList;

   delete m_pObsFiles;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
InitFromFile()

Reads the observation data for each observation point from the input file.
******************************************************************************/
void ObservationGroup::InitFromFile(IroncladString obsFileName)
{
   int i, j; 
   char tmpName[DEF_STR_SZ];
   //int  ID;
   double value;
   double weight;
   char tmpKey[DEF_STR_SZ];
   char group[DEF_STR_SZ];
   int line;
   int col;
   char tmpFile[DEF_STR_SZ];
   char * lineStr;
   char * pTok, tok;
   char tmp1[DEF_STR_SZ];
   FILE * pObsFile;
   bool bQuitOnError;
   double errorVal;

   pObsFile = fopen(obsFileName, "r");

   if(pObsFile == NULL)
   {
      FileOpenFailure("ObservationGroup::InitFromFile", obsFileName);
   }/* end if() */

   //make sure correct tokens are present
   FindToken(pObsFile, "BeginObservations", obsFileName);
   FindToken(pObsFile, "EndObservations", obsFileName);
   rewind(pObsFile);

   //check for optional error-handling instructions
   bQuitOnError = true;
   errorVal = 0.00;
   if(CheckToken(pObsFile, "OnObsError", obsFileName) == true)
   {
      lineStr = GetCurDataLine();
      MyStrLwr(lineStr);
      if(strstr(lineStr, "quit") == NULL)
      {
         sscanf(lineStr, "%s %lf", tmp1, &errorVal);
         bQuitOnError = false;
      }
   }
   rewind(pObsFile);

   //count number of observations
   m_NumObs = 0;
   FindToken(pObsFile, "BeginObservations", obsFileName);
   lineStr = GetNxtDataLine(pObsFile, obsFileName);

   while(strstr(lineStr, "EndObservations") == NULL)
   {   
      m_NumObs++;
      lineStr = GetNxtDataLine(pObsFile, obsFileName);
   }/* end while() */
   rewind(pObsFile);

   //read in each observation
   i = 0;
   NEW_PRINT("Observation *", m_NumObs);
   m_pObsList = new Observation * [m_NumObs];

   if(m_NumObs == 0)
   {
      LogError(ERR_FILE_IO,"No observations specified");
      ExitProgram(1);
   }/* end if() */

   FindToken(pObsFile, "BeginObservations", obsFileName);
   lineStr = GetNxtDataLine(pObsFile, obsFileName);
   while(strstr(lineStr, "EndObservations") == NULL)
   {      
      pTok = lineStr;
      //extract name of observation (no spaces allowed)
      j = ExtractString(pTok, tmpName);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      pTok += j;
      //extract value      
      j = ExtractString(pTok, tmp1);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      value = atof(tmp1);
      pTok += j;
      //extract weight
      j = ExtractString(pTok, tmp1);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      weight = atof(tmp1);
      pTok += j;
      //extract filename (spaces allowed)
      j = ExtractFileName(pTok, tmpFile);
      pTok += j;
      //extract keyword
      j = ExtractString(pTok, tmpKey);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      pTok += j;
      //extract line
      j = ExtractString(pTok, tmp1);
      j = ValidateExtraction(j, 1, 1, "ObservationGroup()");
      line = atoi(tmp1);
      pTok += j;
      //extract column
      j = ExtractString(pTok, tmp1);
      col = atoi(tmp1);
      if(j < 0)
         pTok = &(pTok[strlen(pTok)]);
      else
         pTok += j;
      //extract token (wrapped in ' chars)
      j = ExtractString(pTok, tmp1);
      if((tmp1[0] == 0x27) && (tmp1[2] == 0x27)){ tok = tmp1[1];}
      else{ tok = ' ';}
      if(j < 0)
         pTok = &(pTok[strlen(pTok)]);
      else
         pTok += j;
      if((tok == ' ')  && (tmp1[0] == 0x27)){ pTok+=2;}
      if((tok == '\t') && (tmp1[0] == 0x27)){ pTok+=2;}

      //extract augmented output flag
      bool bAug;
      j = ExtractString(pTok, tmp1);
      if(strcmp(tmp1, "yes") == 0) bAug = true;
      else{ bAug = false;}
      if(j < 0)
         pTok = &(pTok[strlen(pTok)]);
      else
         pTok += j;

      //extract observation group
      strcpy(group, "");
      j = ExtractString(pTok, group);
      if(strcmp(group, "") == 0) strcpy(group, "none");
      if(j < 0)
         pTok = &(pTok[strlen(pTok)]);
      else
         pTok += j;

      //check weight, if it is zero and not part of augmented output, then ignore the observation
      if((fabs(weight) <= NEARLY_ZERO) && (bAug == false))
      {
         sprintf(tmp1, "%s has a weight of zero and has been excluded from the calibration", tmpName);
         LogError(ERR_BAD_WGHT, tmp1);
      }
      else
      {
         NEW_PRINT("Observation", 1);
         m_pObsList[i] = new Observation(tmpName,value,weight,tmpFile,tmpKey,
                                         line,col, tok, bAug, group);

         /*------------------------------------------------------------------
         Create a ValueExtractor class for the given file (if one with that 
         name (i.e. tmpFile) hasn't already been created).
         ------------------------------------------------------------------*/      
         if(m_pObsFiles == NULL)
         { 
            NEW_PRINT("ValueExtractor", 1);
            m_pObsFiles = new ValueExtractor(tmpFile, bQuitOnError, errorVal);
         }
         else{ m_pObsFiles->Insert(tmpFile);}

         i++;      
      }
      lineStr = GetNxtDataLine(pObsFile, obsFileName);
   }/* end while() */

   if(i < m_NumObs)
   {
      for(j = i; j < m_NumObs; j++) m_pObsList[j] = NULL;
      m_NumObs = i;
   }
      
   fclose(pObsFile);

   m_NumGroups = GetNumGroups();
   /*
   for(i = 0; i < m_NumGroups; i++)
   {
      printf("Group %d Name = %s\n", i, GetGroup(i));
   }*/
} /* end InitFromFile() */

/******************************************************************************
ExcludeObs()

Remove the given observation from the active list.
******************************************************************************/
void ObservationGroup::ExcludeObs(UnchangeableString obs)
{
   int i, j;

   for(i = 0; i < m_NumObs; i++)
   {
      if(strcmp(obs, m_pObsList[i]->GetName()) == 0) break;
   }

   if(i == m_NumObs) return; //no match

   delete m_pObsList[i];

   for(j = i; j < (m_NumObs-1); j++)
   {
      m_pObsList[j] = m_pObsList[j+1];
   }/* end for() */
   m_pObsList[j] = NULL;
   m_NumObs--;
}/* end ExcludeObs() */

/******************************************************************************
Write()

Writes user-specified simulated output to the pFile argument.
******************************************************************************/
void ObservationGroup::Write(FILE * pFile, int type, double * F)
{
   int i;
   
   //emit MO data
   if(AlgIsMultiObjective() == true)
   {
      for(i = 0; i < m_NumGroups; i++)
      {
         if(type == WRITE_BNR)
         {
            fprintf(pFile, "WSSE(%-6s)  ", GetGroup(i));
         }
         else
         {
            fprintf(pFile, "%E  ", F[i]);
         }
      }/* end for() */
   }/* end if() */

   for(i = 0; i < m_NumObs; i++)
   {
      if(m_pObsList[i]->IsAugmented())
      {
         m_pObsList[i]->WriteSim(pFile, type);
      }
   }
} /* end write() */
