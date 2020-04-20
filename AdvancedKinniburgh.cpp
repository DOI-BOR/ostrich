/******************************************************************************
File     : AdvancedKinniburgh.cpp
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Orear isotherm model. Provides a text interface for the set of Isotherms that 
make up the IsoFit and Isotherm programs. Unlike IsoFit.cpp, AdvancedKinniburgh.cpp only
uses measured 'C' values and does not require measured 'q' values. The solution 
method is based on:

David G. Kinniburgh. 1986. General Purpose Adsorption Isotherms. Environmental 
Science and Technology, vol. 20, no. 9, pg. 895-904.

Here is an explanation of the Kinniburgh approach, with much text borrowed from 
Alan J. Rabideau (personal communication):

Kinniburgh proposes that we should estimate isotherm parameters by minimizing 
the weighted sum of squares of the aqueous concentrations, which are predicted 
as follows:
   C = Ct/(1+X) - [(S/V)/(1+X)] * q(C)
where Ct = the initial aqueous concentration in the vial
   S = solid/liquid ratio in vial (SolidLiquidRatio)
   V = volume of liquid in vial (TotalVolume)
   q = the isotherm expression (q calculated from C)
   X = proportional loss term
This formulation should resolve the measurement error problem, since errors
in Ct, S, and V can be assumed small (much smaller than C). See Eqs. 6-8 in 
Kinniburgh (1986).

Since 'C' appear on both sides of the equation, the Kinniburgh approach is 
non-linear. Also note that the Kinniburgh formulation avoids measured q values.  
The reasoning is that since q is  usually computed from C, rather than measured, 
it is only an indirect measure of the "real" dependent variable.  The manner of 
estimating q errors by error propagation necessitate that they are larger than the 
errors in C, supporting the common the justification for casting the residuals in 
terms of C; that is, the errors in q are assumed sufficiently large to overwhelm 
the errors in C, which allows traditional regression methods to be used.

If we actually measured q (eg, by separating the solid phase and digesting it), 
we would then have the "measurement errors in both variables" situation. In such a 
case, user should consider using either the McCammon or Orear solution methods, as 
these explicitly deal with error in both variables.

The Kinniburgh approach is nonlinear and requires the following steps: 
      1) Assign initial isotherm parameters
      2) For each data point, Calculate C that minimizes the following 
         nonlinear equation:
         minimize |C - [Ct/(1+X) - [(S/V)/(1+X)] * q(C)]|^2
      3) calculate WSSE as a function of the measured and calculated C's
      4) update isotherm parameters
      5) return to step 2

The Kinniburgh() routine is only responsible for Step 2. We let the Ostrich search
algorithm tackle steps 1 and 3-5 as part of the Isotherm parameter estimation.

Version History
03-10-10    lsm   Created, based on Kinniburgh.cpp
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "AdvancedKinniburghSolver.h"
#include "Isotherms.h"

#include "IsoParse.h"
#include "Exception.h"

int AdvancedKinniburgh(void)
{  
   AdvancedKinniburghSolver * pKini;
   int size;
   char * pLine;
   char pVar[DEF_STR_SZ];
   char pType[DEF_STR_SZ];
   char * pStr;
   char * pTmp;
   IsothermABC * pIso;

   //printf("Advanced Kinniburgh Program\n");

   size = ISO_GetFileSize(ISO_IN_FILE);
   if(size <= 0)
   {
      LogError(ERR_FILE_IO, "AdvancedKinniburgh() : empty or nonexistant input file");
      ExitProgram(1);
   }
   NEW_PRINT("char", size+1);
   pStr = new char[size+1];

   MEM_CHECK(pStr);
   ISO_FileToStr(ISO_IN_FILE, pStr, size);

   pTmp = strstr(pStr, "IsothermType");
   if(pTmp == NULL)
   {
      LogError(ERR_BAD_ARGS, "AdvancedKinniburgh() : Unspecified isotherm type");
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
      LogError(ERR_BAD_ARGS, "AdvancedKinniburgh() : Unknown isotherm type, valid types are:");
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
      LogError(ERR_FILE_IO, "AdvancedKinniburgh() : could not initialize Isotherm");
      delete [] pStr;
      delete pIso;
      ExitProgram(1);
   }

   //extract X parameter
   double X;
   pTmp = strstr(pStr, "KinniburghLossTerm");
   if(pTmp == NULL) X = 1.00;
   else sscanf(pTmp, "KinniburghLossTerm %lf", &X);

   pKini = new AdvancedKinniburghSolver(pIso, X);
   MEM_CHECK(pKini);
   //printf("Initializing Kinniburgh Solver\n");
   if(pKini->Initialize(pStr) == false)
   {
      LogError(ERR_FILE_IO, "AdvancedKinniburgh() : could not initialize solver");
      delete [] pStr;
      delete pIso;
      delete pKini;
      ExitProgram(1);
   }
   
   //printf("Kinniburgh() Solving non-linear expressions for Ci\n");
   pKini->Compute();

   //printf("Done. Output stored in: |%s|\n", ISO_OUT_FILE);
   delete [] pStr;
   delete pIso;
   delete pKini;
   return (0);
} /* end AdvancedKinniburgh() */
