/******************************************************************************
File     : McCammon.cpp
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

McCammon isotherm model. Provides a text interface for the set of Isotherms that 
make up the IsoFit and Isotherm programs. Unlike IsoFit.cpp, McCammon.cpp handles
errors-in-C. Solution method is based on:

RB McCammon. 1973. Nonlinear Regression for Dependent Variables. Mathematical
Geology, vol. 5, no. 4, pg. 365-375.

The cornerstone of the McCammon method is the following non-linear equation (see
McCammon 1973, equation 10):

   (Cobs - Cest)+ dq*[wq*wq/wc*wc]*[qobs - q(Cest)] = 0
where:
   Cest is the simulated aqueous concentration
   Cobs is the measured aqueous concentration
   qobs is the measured sorbed concentration
   q(Cobs) is the simulated sorbed concentration (i.e. the isotherm expression)
   wc is the aqueous observation weight
   wq is the sorbed observation weight
   dq is the derivative of q() [i.e. dq/dc], evaluated at Cest

Since the equation in non-linear with respect to Cest we must iterate to find 
a solution. 
      1) Assign initial isotherm parameters
      2) For each data point, Calculate Cest that minimizes the following 
         nonlinear equation:
         minimize [F(Cest)]^2, where F() is the eqn. given above
      3) calculate WSSE as a function of the measured and simulated C's
      4) update isotherm parameters
      5) return to step 2

The McCammon() routine is only responsible for Step 2. We let the Ostrich search
algorithm tackle steps 1 and 3-5 as part of the Isotherm parameter estimation.

Version History
07-28-07    lsm   Created, based on Orear.cpp
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "McCammonSolver.h"
#include "Isotherms.h"

#include "IsoParse.h"
#include "Exception.h"

IsothermABC * gMcIso = NULL;
McCammonSolver * gMcCam = NULL;

void DisklessMcCammon(ParameterGroup * pgroup, ObservationGroup * ogroup)
{
   if((gMcIso == NULL) || (gMcCam == NULL)){ return; }
   if((pgroup == NULL) && (ogroup == NULL)) { delete gMcCam; delete gMcIso; return; }
   gMcIso->Initialize(pgroup);
   gMcCam->Compute(ogroup);
}/* end DisklessIsotherm() */

int McCammon(bool bSave)
{  
   McCammonSolver * pMcCam;
   int size;
   char * pLine;
   char pVar[DEF_STR_SZ];
   char pType[DEF_STR_SZ];
   char * pStr;
   char * pTmp;
   IsothermABC * pIso;

   // printf("McCammon Program\n");

   size = ISO_GetFileSize(ISO_IN_FILE);
   if(size <= 0)
   {
      LogError(ERR_FILE_IO, "McCammon() : empty or nonexistant input file");
      ExitProgram(1);
   }
   NEW_PRINT("char", size+1);
   pStr = new char[size+1];

   MEM_CHECK(pStr);
   ISO_FileToStr(ISO_IN_FILE, pStr, size);

   pTmp = strstr(pStr, "IsothermType");
   if(pTmp == NULL)
   {
      LogError(ERR_BAD_ARGS, "McCammon() : Unspecified isotherm type");
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
      LogError(ERR_BAD_ARGS, "McCammon() : Unknown isotherm type, valid types are:");
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
      LogError(ERR_FILE_IO, "McCammon() : could not initialize Isotherm");
      delete [] pStr;
      delete pIso;
      ExitProgram(1);
   }

   pMcCam = new McCammonSolver(pIso);
   MEM_CHECK(pMcCam);
   //printf("Initializing McCammon Solver\n");
   if(pMcCam->Initialize(pStr) == false)
   {
      LogError(ERR_FILE_IO, "McCammon() : could not initialize solver");
      delete [] pStr;
      delete pIso;
      delete pMcCam;
      ExitProgram(1);
   }
   
   //printf("McCammon() Solving non-linear expressions\n");
   pMcCam->Compute();

   //printf("Done. Output stored in: |%s|\n", ISO_OUT_FILE);
   delete [] pStr;

   if(bSave == true)
   {
      gMcIso = pIso;
      gMcCam = pMcCam;
   }
   else
   {
      delete pIso;
      delete pMcCam;
   }

   return (0);
} /* end McCammon() */
