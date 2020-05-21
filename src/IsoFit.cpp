/******************************************************************************
File     : IsoFit.cpp
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Isoterhm model. Provides a text interface for the set of Isotherms that 
make up the IsoFit and Isotherm programs.

Version History
03-16-04    lsm   Created
03-24-04    lsm   Removed printf() statements
08-17-04    lsm   Added printing of RAM allocations
03-04-05    lsm   Added Polanyi-Partition isotherm
06-09-05    lsm   Added Langmuir-partition isotherm
01-01-07    lsm   Added additional isotherms, the current list of options are:
                     1.  BET Isotherm
                     2.  Freundlich Isotherm
                     3.  Freundlich-Partition Isotherm
                     4.  Linear Isotherm
                     5.  Langmuir Isotherm
                     6.  Generalized Langmuir-Freundlich Isotherm
                     7.  Langmuir-Partition Isotherm
                     8.  Polanyi Isotherm
                     9.  Polanyi-Partition Isotherm
                     10. Toth Isotherm
06-22-07    lsm      11. Dual-Langmuir Isotherm
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Isotherms.h"

#include "IsoParse.h"
#include "Exception.h"

IsothermABC * gIso = NULL;

void DisklessIsotherm(ParameterGroup * pgroup, ObservationGroup * ogroup)
{
   if(gIso == NULL){ return; }
   if((pgroup == NULL) && (ogroup == NULL)) { delete gIso; return; }
   gIso->Initialize(pgroup);
   gIso->Compute(ogroup);
}/* end DisklessIsotherm() */

int Isotherm(bool bSave)
{  
   int size;
   char * pLine;
   char pVar[DEF_STR_SZ];
   char pType[DEF_STR_SZ];
   char * pStr;
   char * pTmp;
   IsothermABC * pIso;

   //printf("Isotherm Program\n");

   size = ISO_GetFileSize(ISO_IN_FILE);
   if(size <= 0)
   {
      LogError(ERR_FILE_IO, "Isotherm() : empty or nonexistant input file");
      ExitProgram(1);
   }
   NEW_PRINT("char", size+1);
   pStr = new char[size+1];
   MEM_CHECK(pStr);

   ISO_FileToStr(ISO_IN_FILE, pStr, size);

   pTmp = strstr(pStr, "IsothermType");
   if(pTmp == NULL)
   {
      LogError(ERR_BAD_ARGS, "Isotherm() : Unspecified isotherm type");
      delete [] pStr;
      ExitProgram(1);
   }
   ISO_GetLine(pTmp, &pLine);
   sscanf(pLine, "%s %s", pVar, pType);

   if(strcmp(pType, "LinearIsotherm") == 0)
   { 
      NEW_PRINT("LinearIsotherm", 1); 
      pIso = new LinearIsotherm();
   }
   else if(strcmp(pType, "LangmuirIsotherm") == 0)
   {
      NEW_PRINT("LangmuirIsotherm", 1); 
      pIso = new LangmuirIsotherm();
   }
   else if(strcmp(pType, "DualLangmuirIsotherm") == 0)
   {
      NEW_PRINT("DualLangmuirIsotherm", 1); 
      pIso = new DualLangmuirIsotherm();
   }
   else if(strcmp(pType, "FreundlichIsotherm") == 0)
   { 
      NEW_PRINT("FreundlichIsotherm", 1);
      pIso = new FreundlichIsotherm();
   }
   else if(strcmp(pType, "Polanyi-PartitionIsotherm") == 0)
   { 
      NEW_PRINT("PolanyiPartitionIsotherm", 1);
      pIso = new PolanyiPartitionIsotherm();
   }
   else if(strcmp(pType, "Langmuir-PartitionIsotherm") == 0)
   { 
      NEW_PRINT("LangmuirPartitionIsotherm", 1);
      pIso = new LangmuirPartitionIsotherm();
   }
   else if(strcmp(pType, "BET_Isotherm") == 0)
   { 
      NEW_PRINT("BET_Isotherm", 1);
      pIso = new BET_Isotherm();
   }
   else if(strcmp(pType, "TothIsotherm") == 0)
   { 
      NEW_PRINT("TothIsotherm", 1);
      pIso = new TothIsotherm();
   }
   else if(strcmp(pType, "Langmuir-FreundlichIsotherm") == 0)
   { 
      NEW_PRINT("LangmuirFreundlichIsotherm", 1);
      pIso = new LangmuirFreundlichIsotherm();
   }
   else if(strcmp(pType, "PolanyiIsotherm") == 0)
   { 
      NEW_PRINT("PolanyiIsotherm", 1);
      pIso = new PolanyiIsotherm();
   }
   else if(strcmp(pType, "Freundlich-PartitionIsotherm") == 0)
   { 
      NEW_PRINT("FreundlichPartitionIsotherm", 1);
      pIso = new FreundlichPartitionIsotherm();
   }
   else if(strcmp(pType, "OrearIsotherm") == 0)
   {
      NEW_PRINT("OrearIsotherm", 1); 
      pIso = new OrearIsotherm();
   }
   else if(strcmp(pType, "McCammonIsotherm") == 0)
   {
      NEW_PRINT("McCammonIsotherm", 1); 
      pIso = new McCammonIsotherm();
   }
   else
   {
      LogError(ERR_BAD_ARGS, "Isotherm() : Unknown isotherm type, valid types are:");
      LogError(ERR_CONTINUE, "**********************************");
      LogError(ERR_CONTINUE, "   BET_Isotherm");
      LogError(ERR_CONTINUE, "   FreundlichIsotherm");
      LogError(ERR_CONTINUE, "   Freundlich-PartitionIsotherm");
      LogError(ERR_CONTINUE, "   LinearIsotherm");
      LogError(ERR_CONTINUE, "   LangmuirIsotherm");
      LogError(ERR_CONTINUE, "   DualLangmuirIsotherm");
      LogError(ERR_CONTINUE, "   Langmuir-FreundlichIsotherm");
      LogError(ERR_CONTINUE, "   Langmuir-PartitionIsotherm");
      LogError(ERR_CONTINUE, "   PolanyiIsotherm");
      LogError(ERR_CONTINUE, "   Polanyi-PartitionIsotherm");
      LogError(ERR_CONTINUE, "   TothIsotherm");
      LogError(ERR_CONTINUE, "**********************************");

      delete [] pStr;
      ExitProgram(1);
   }
   MEM_CHECK(pIso);

   //printf("Initializing %s\n", pType);
   if(pIso->Initialize(pStr) == false)
   {
      LogError(ERR_FILE_IO, "Isotherm() : could not initialize Isotherm");
      delete [] pStr;
      delete pIso;
      ExitProgram(1);
   }
   //printf("Computing sorbed concentrations (qi)\n");
   pIso->Compute();
   //printf("Done. Output stored in: |%s|\n", ISO_OUT_FILE);
   delete [] pStr;

   if(bSave == true) gIso = pIso;
   else delete pIso;

   return (0);
} /* end Isotherm() */
