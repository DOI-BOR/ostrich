/******************************************************************************
File     : IsoParse.cpp
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Parsing routines for the IsoTherm model, an internal model that computes 
sorption istherms.

Version History
03-16-04    lsm   Created
03-24-04    lsm   Integrated IsoFit, Isotherm() and Ostrich programs
08-17-04    lsm   Added printing of RAM allocations
03-04-05    lsm   Added Polanyi-Partition isotherm
10-19-05    lsm   Added a Ranges section (to specifiy lower, upper limits and
                  transformations). Added support for user-specified Swarm Size 
                  and Max Generations. Added Langmuir-Partition isotherm.
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
01-01-07    lsm   Added a "fit-all" option that will apply IsoFit to all 
                  isotherms. To Select this option, the user should include the
                  following line in the IsoFitIn.txt input file:
                     IsothermType     AllIsotherms
01-01-07    lsm   Added additional weighting schemes, the current list of options are:
                     1. Uniform <weight>
                     2. SorbedRelative  <relative error>
                     3. AqueousRelative <relative error>  <conversion factor>
                     4. IndividualStdDevs
01-01-07    lsm   Added additional stat calculations, currently computed stats are:
                    StdDev and StdErr : basic goodness of fir measures
                    CorrCoeff : parameter correlation
                    Beale and Linssen : measures of isotherm linearity
                    CooksD and DFBETAS : measures of observation influence
                    Confidence : linear CI for parameters
                    NormPlot : R-squared for normally distributed residuals
                    Sensitivity : parameter sensitivity
                    Matrices : the Jacobian/Sensitivity matrix
                    RunsTest and Durbin-Watson : for testing serial autocorrelation of residuals
                    MMRI : for comparing multiple isotherms
05-25-06    lsm   Added additional isotherms:
                     11. Dual Langmuir Isotherm
07-29-07    lsm   Added support for Kinniburgh solution method
08-28-07    lsm   Added support for holding (or not) insensitive observations and paramaters
03-10-10    lsm   Added support for Advanced Kinniburgh solution method
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ParameterGroup.h"
#include "ObservationGroup.h"

#include "IsoParse.h"
#include "Exception.h"
#include "Utility.h"
#include "StatUtility.h"

#define DEFAULT_STEP (0.001)
#define NUM_ISOTHERMS (11)

/*
A listing of isotherms for the "fit all" options.
*/
IroncladString gIsoNames[NUM_ISOTHERMS] = 
{
   "BET_",
   "DualLangmuir",
   "Freundlich",
   "Freundlich-Partition",
   "Linear",
   "Langmuir",
   "Langmuir-Freundlich",
   "Langmuir-Partition",
   "Polanyi",
   "Polanyi-Partition",
   "Toth"
};

char gBeginRanges[DEF_STR_SZ] = "BeginRanges";
char gEndRanges[DEF_STR_SZ]   = "EndRanges";
char gISO_PSOOUT_FILE[DEF_STR_SZ];
char gISO_PSOMDL_FILE[DEF_STR_SZ];

/* 
Global storage for line of input.
*/
char * gIsoLine = NULL;
int gIsoLineSize = 0;

/******************************************************************************
main()

Entry point for the IsoFit application, which uses the Isotherm internal model 
and Ostrich routines to calibrate isotherm data.
******************************************************************************/
#ifdef ISOFIT_BUILD
//#warning error
int main(int argc, StringType argv[])
{
   NEW_PRINT("IsoParamList", 1);
   IsoParamList * pList;
   IsoGlobStruct Glob;

   Glob.conc = NULL;
   Glob.sorb = NULL;
   Glob.wsorb = NULL;
   Glob.wconc = NULL;
   Glob.expA = NULL;
   Glob.expB = NULL;
   Glob.expD = NULL;

   char tmpStr[DEF_STR_SZ];
   bool bDone = false;
   int i = 0;

   InitErrors();

   ISO_ReadIsoFitFile(&Glob);
   if(Glob.numObs == 0)
   {
      LogError(ERR_FILE_IO, "Error reading IsoFit input file");
      ExitProgram(1);
   }

   while(bDone == false)
   {
      strcpy(gISO_PSOOUT_FILE, ISO_PSOOUT_FILE);
      strcpy(gISO_PSOMDL_FILE, ISO_PSOMDL_FILE);

      if(Glob.bFitAll == true)
      { 
         strcpy(Glob.isoStr, gIsoNames[i]); 
         strcat(Glob.isoStr, "Isotherm");

         strcpy(gBeginRanges, "Begin");
         strcat(gBeginRanges, gIsoNames[i]);
         strcat(gBeginRanges, "Ranges");

         strcpy(gEndRanges, "End");
         strcat(gEndRanges, gIsoNames[i]);
         strcat(gEndRanges, "Ranges");

         strcat(gISO_PSOOUT_FILE, "_");
         strcat(gISO_PSOMDL_FILE, "_");
         strcat(gISO_PSOOUT_FILE, gIsoNames[i]);
         strcat(gISO_PSOMDL_FILE, gIsoNames[i]);
         strcat(gISO_PSOOUT_FILE, ".txt");
         strcat(gISO_PSOMDL_FILE, ".txt");
      }
      else
      {
         strcat(gISO_PSOOUT_FILE, ".txt");
         strcat(gISO_PSOMDL_FILE, ".txt");

         bDone = true;
      }

      sprintf(tmpStr, "OstOutput0_%s.txt", gIsoNames[i]);
      if(MY_ACCESS(tmpStr, 0) != -1)
      {
        printf("File %s exists, skipping Isotherm fitting\n", tmpStr);
        i++;
        if(i >= NUM_ISOTHERMS){ bDone = true;}
      }
      else
      {
         ISO_CreateTemplateFile(&Glob);

         pList = new IsoParamList;

         ISO_CreateParamList(pList, &Glob);

         ISO_CreateOstrichFile("ParticleSwarm", pList, false, &Glob);

         /* perform user-defined Ostrich pre-processing */
         if(MY_ACCESS("OstrichPreProcessor.bat", 0) != -1)
         {
            system("OstrichPreProcessor.bat");
         }

         Ostrich(argc, argv);
         ISO_RefreshParamList(pList);
         remove(gISO_PSOOUT_FILE);
         remove(gISO_PSOMDL_FILE);
         rename(ISO_OSTOUT_FILE, gISO_PSOOUT_FILE);
         rename(ISO_OSTMDL_FILE, gISO_PSOMDL_FILE);

         ISO_CreateOstrichFile("Powell", pList, false, &Glob);
         Ostrich(argc, argv);
         ISO_RefreshParamList(pList);

         ISO_CreateOstrichFile("Levenberg-Marquardt", pList, true, &Glob);
         Ostrich(argc, argv);

         ISO_DestroyIsoParamList(pList);
    
         if(Glob.debug == 0)
         {
            remove(gISO_PSOMDL_FILE);
            remove(gISO_PSOOUT_FILE);
         }/* end if() */

         //When in "fit all" mode, must rename files to prevent them from being overwritten
         if(Glob.bFitAll == true)
         {
            sprintf(tmpStr, "OstOutput0_%s.txt", gIsoNames[i]);
            remove(tmpStr);
            rename("OstOutput0.txt", tmpStr);

            sprintf(tmpStr, "OstModel0_%s.txt", gIsoNames[i]);
            remove(tmpStr);
            rename("OstModel0.txt", tmpStr);

            sprintf(tmpStr, "OstJacobian0_%s.txt", gIsoNames[i]);
            remove(tmpStr);
            rename("OstJacobian0.txt", tmpStr);

            ReportErrors();
            sprintf(tmpStr, "OstErrors0_%s.txt", gIsoNames[i]);
            remove(tmpStr);
            rename("OstErrors0.txt", tmpStr);

            sprintf(tmpStr, "OstStatus0_%s.txt", gIsoNames[i]);
            remove(tmpStr);
            rename("OstStatus0.txt", tmpStr);

            sprintf(tmpStr, "%s_%s.txt", GetOstExeOut(), gIsoNames[i]);
            remove(tmpStr);
            rename(GetOstExeOut(), tmpStr);

            if(Glob.debug == 1)
            {
               sprintf(tmpStr, "IsothermIn_%s.tpl", gIsoNames[i]);
               remove(tmpStr);
               rename("IsothermIn.tpl", tmpStr);

               sprintf(tmpStr, "IsothermIn_%s.txt", gIsoNames[i]);
               remove(tmpStr);
               rename("IsothermIn.txt", tmpStr);

               sprintf(tmpStr, "IsothermOut_%s.txt", gIsoNames[i]);
               remove(tmpStr);
               rename("IsothermOut.txt", tmpStr);

               sprintf(tmpStr, "OstIn_%s.txt", gIsoNames[i]);
               remove(tmpStr);
               rename("OstIn.txt", tmpStr);
            }/* end if() */
            i++;
            if(i >= NUM_ISOTHERMS){ bDone = true;}
         }/* end if(bFitAll == true) */
      }/* end else() */
   }/* end while(bDone == false) */
   delete [] Glob.conc;
   delete [] Glob.sorb;
   delete [] Glob.wsorb;
   delete [] Glob.wconc;
   delete [] Glob.expA;
   delete [] Glob.expB;
   delete [] Glob.expD;
   
   if(Glob.debug == 0)
   {
      remove("IsothermIn.tpl");
      remove("IsothermIn.txt");
      remove("IsothermOut.txt");
      remove("OstIn.txt");
   }/* end if() */
   ExitProgram(0);
   return 0;
}/* end main() */
#endif

/******************************************************************************
ISO_GetSwarmParams()

Read in the max. generations and population size to use in the particle swarm
optimization procedure.
******************************************************************************/
void ISO_GetSwarmParams(char * pStr, IsoGlobStruct * pArgs)
{
   char * pTmp = NULL;   
   char * pLine;
   bool bPop = false, bGen = false;

   //check tokens
   pTmp = pStr;

   /*----------------------------------------------
   Read In Parameters
   ----------------------------------------------*/
   pArgs->popSize = 0;
   pArgs->maxGens = 0;
   do{ 
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strncmp(pLine, "PopSize", 7) == 0)
      {
         pArgs->popSize = atoi(&(pLine[7]));
         bPop = true;
      }
      if(strncmp(pLine, "MaxGens", 7) == 0)
      {
         pArgs->maxGens = atoi(&(pLine[7]));
         bGen = true;
      }
      if((bPop == true) && (bGen == true)) break;
   }while(*pTmp != (char)NULL);
}/* end ISO_GetSwarmParams() */

/******************************************************************************
ISO_GetMethod()

Read in the solution method to employ for isotherm fitting.
******************************************************************************/
IsoFitSolverType ISO_GetMethod(char * pStr)
{
   char * pTmp = NULL;   
   char * pLine, pTok[DEF_STR_SZ], pChoice[DEF_STR_SZ];
   IsoFitSolverType method;

   //check tokens
   pTmp = pStr;

   /*----------------------------------------------
   Read in user setting, if one exists
   ----------------------------------------------*/
   method = ISOTHERM_METHOD;
   do
   { 
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strncmp(pLine, "SolutionMethod", 14) == 0)
      {
         sscanf(pLine, "%s %s", pTok, pChoice);
         if(strcmp(pChoice, "Standard") == 0)
         {
            method = ISOTHERM_METHOD;
         }
         else if (strcmp(pChoice, "Orear") == 0)
         {
            method = OREAR_METHOD;
         }
         else if (strcmp(pChoice, "McCammon") == 0)
         {
            method = MCCAMMON_METHOD;
         }
         else if (strcmp(pChoice, "Kinniburgh") == 0)
         {
            method = KINNIBURGH_METHOD;
         }
         else if (strcmp(pChoice, "AdvancedKinniburgh") == 0)
         {
            method = ADV_KINNIBURGH_METHOD;
         }
         else if (strcmp(pChoice, "TotalError") == 0)
         {
            method = TOTAL_ERROR_METHOD;
         }
         else
         {
            LogError(ERR_BAD_ARGS, "ISO_GetMethod() : Unknown method, valid methods are:");
            LogError(ERR_CONTINUE, "**********************************");
            LogError(ERR_CONTINUE, "   Standard");
            LogError(ERR_CONTINUE, "   Orear");
            LogError(ERR_CONTINUE, "   McCammon");
            LogError(ERR_CONTINUE, "   Kinniburgh");
            LogError(ERR_CONTINUE, "   TotalError");
            LogError(ERR_CONTINUE, "   AdvancedKinniburgh");
            LogError(ERR_CONTINUE, "**********************************");

            delete [] pStr;
            ExitProgram(1);
         }
         break;
      }
   }while(*pTmp != (char)NULL);

   return method;
}/* end ISO_GetMethod() */

/******************************************************************************
ISO_ReadIsoFitFile()

Read in the name of the IsoTherm, the list of aqueous and sorbed concentrations, 
and the weighting scheme.

Returns number of observations (size of conc, sorb and weight arrays)
******************************************************************************/
void ISO_ReadIsoFitFile(IsoGlobStruct * pArgs)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   char pType[DEF_STR_SZ];
   double val;
   double convFactor;
   int size;
   char * pStr;
   double sd, xsd;
   FILE * pFile;

   pArgs->numObs = 0; //signals an error

   size = ISO_GetFileSize("IsoFitIn.txt");
   if(size <= 0)
   {
      LogError(ERR_FILE_IO, "ISO_ReadIsoFitFile() : empty or nonexistant input file");
      ExitProgram(1);
   }
   NEW_PRINT("char", size+1);
   pStr = new char[size+1];
   MEM_CHECK(pStr);
   ISO_FileToStr("IsoFitIn.txt", pStr, size);

   //check tokens
   msg[0] = (char)NULL;
   if((pTmp=strstr(pStr, "BeginLabData")) == NULL){ strcat(msg, "BeginLabData\n");}
   if(strstr(pTmp, "EndLabData")          == NULL){ strcat(msg, "EndLabData\n");}
   if(strstr(pStr, "WeightingScheme") == NULL){ strcat(msg, "WeightingScheme\n");}
   if(strstr(pStr, "IsothermType") == NULL){ strcat(msg, "IsothermType\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      delete [] pStr;
      ExitProgram(1);
   }

   /*-------------------------------------------------------------------------------------
   Determine Isotherm type
   -------------------------------------------------------------------------------------*/
   pTmp = strstr(pStr, "IsothermType");   
   ISO_GetLine(pTmp, &pLine);
   sscanf(pLine, "%s %s", pVar, pType);
   if((strcmp(pType, "AllIsotherms")                 == 0) ||
      (strcmp(pType, "BET_Isotherm")                 == 0) ||
      (strcmp(pType, "DualLangmuirIsotherm")         == 0) ||
      (strcmp(pType, "FreundlichIsotherm")           == 0) ||
      (strcmp(pType, "Freundlich-PartitionIsotherm") == 0) ||
      (strcmp(pType, "LinearIsotherm")               == 0) ||
      (strcmp(pType, "LangmuirIsotherm")             == 0) ||
      (strcmp(pType, "Langmuir-FreundlichIsotherm")  == 0) ||
      (strcmp(pType, "Langmuir-PartitionIsotherm")   == 0) ||
      (strcmp(pType, "PolanyiIsotherm")              == 0) ||
      (strcmp(pType, "Polanyi-PartitionIsotherm")    == 0) ||
      (strcmp(pType, "TothIsotherm")                 == 0) ||
      (strcmp(pType, "OrearIsotherm")                == 0) ||
      (strcmp(pType, "McCammonIsotherm")             == 0))
   { strcpy(pArgs->isoStr, pType);}
   else
   {
      LogError(ERR_BAD_ARGS, "ISO_ReadIsoFitFile() : Unknown isotherm type, valid types are:");
      LogError(ERR_CONTINUE, "**********************************");
      LogError(ERR_CONTINUE, "   AllIsotherms");
      LogError(ERR_CONTINUE, "   BET_Isotherm");
      LogError(ERR_CONTINUE, "   DualLangmuirIsotherm");
      LogError(ERR_CONTINUE, "   FreundlichIsotherm");
      LogError(ERR_CONTINUE, "   Freundlich-PartitionIsotherm");
      LogError(ERR_CONTINUE, "   LinearIsotherm");
      LogError(ERR_CONTINUE, "   LangmuirIsotherm");
      LogError(ERR_CONTINUE, "   Langmuir-FreundlichIsotherm");
      LogError(ERR_CONTINUE, "   Langmuir-PartitionIsotherm");
      LogError(ERR_CONTINUE, "   PolanyiIsotherm");
      LogError(ERR_CONTINUE, "   Polanyi-PartitionIsotherm");
      LogError(ERR_CONTINUE, "   TothIsotherm");
      LogError(ERR_CONTINUE, "**********************************");

      delete [] pStr;
      ExitProgram(1);
   }

   /* Polanyi-Partition Isotherm Requires solubility value */
   strcpy(pArgs->solStr, "n/a");
   if((strcmp(pArgs->isoStr, "Polanyi-PartitionIsotherm") == 0) ||
      (strcmp(pArgs->isoStr, "PolanyiIsotherm") == 0) ||
      (strcmp(pArgs->isoStr, "BET_Isotherm") == 0)    ||
      (strcmp(pArgs->isoStr, "AllIsotherms") == 0))
   {
      pTmp = strstr(pStr, "Solubility");
      if(pTmp == NULL)
      { 
         printf("The following token is missing: Solubility \n");
         delete [] pStr;
      }
      else
      {
         pTmp = ISO_GetLine(pTmp, &pLine);
         sscanf(pLine, "%s %s", pVar, pArgs->solStr);
      }
   }

   /* default to no lumping of Q0 and b terms, then override if requested */
   strcpy(pArgs->lumpStr, "no");
   pArgs->bLumpedQ0 = false;
   pTmp = strstr(pStr, "LumpedQ0*b");
   if(pTmp != NULL)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %s", pVar, pArgs->lumpStr);
      MyTrim(pArgs->lumpStr);
      MyStrLwr(pArgs->lumpStr);
   }

   /*-------------------------------------------------------------------------
   Determine weighting scheme
   --------------------------------------------------------------------------*/
   pTmp = strstr(pStr, "WeightingScheme");
   ISO_GetLine(pTmp, &pLine);
   convFactor = 1.00e00;
   sscanf(pLine, "%s %s %lf %lf", pVar, pType, &val, &convFactor);

   //read in solution method
   pArgs->method = ISO_GetMethod(pStr);

   /*--------------------------------------------------------------------------
   Make sure weighting scheme is appropriate for selected solution method. In 
   particular, the "Relative" weighting schemes don't make sense if the Orear, 
   McCammon or TotalError solution method is chosen. Also, the sorbed-relative 
   weighting scheme doesn't make sense if the Kinniburgh or AdvancedKinniburgh
   solution method is chosen.
   --------------------------------------------------------------------------*/
   if((pArgs->method == OREAR_METHOD) || (pArgs->method == MCCAMMON_METHOD) ||
      (pArgs->method == TOTAL_ERROR_METHOD))
   {
      if ((strcmp(pType, "SorbedRelative") == 0) || 
          (strcmp(pType, "AqueousRelative") == 0))
      {
         LogError(ERR_BAD_ARGS, "ISO_ReadIsoFitFile() : Invalid weighting scheme");
         LogError(ERR_CONTINUE, "For Orear/McCammon/TotalError solution methods,");
         LogError(ERR_CONTINUE, "valid schemes are:");
         LogError(ERR_CONTINUE, "*****************************************************");
         LogError(ERR_CONTINUE, "Uniform <weight>");
         LogError(ERR_CONTINUE, "IndividualStdDevs");
         LogError(ERR_CONTINUE, "*****************************************************");
         delete [] pStr;
         ExitProgram(1);
      }
   }
   else if((pArgs->method == KINNIBURGH_METHOD) || (pArgs->method == ADV_KINNIBURGH_METHOD))
   {
      if (strcmp(pType, "SorbedRelative") == 0)
      {
         LogError(ERR_BAD_ARGS, "ISO_ReadIsoFitFile() : Invalid weighting scheme");
         LogError(ERR_CONTINUE, "For the Kinniburgh and AdvancedKinniburgh solution methods,");
         LogError(ERR_CONTINUE, "valid schemes are:");
         LogError(ERR_CONTINUE, "*****************************************************");
         LogError(ERR_CONTINUE, "Uniform <weight>");
         LogError(ERR_CONTINUE, "AqueousRelative <relative error>  <conversion factor>");
         LogError(ERR_CONTINUE, "IndividualStdDevs");
         LogError(ERR_CONTINUE, "*****************************************************");
         delete [] pStr;
         ExitProgram(1);
      }
   }

   /*-------------------------------------------------------------------------------------
   Read in the lab data section
   -------------------------------------------------------------------------------------*/
   pArgs->numObs = 0;
   pTmp = strstr(pStr,"BeginLabData");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strstr(pLine, "EndLabData") == NULL)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      pArgs->numObs++;
   }
   pArgs->numObs--;

   NEW_PRINT("double", pArgs->numObs);
   pArgs->expA = new double[pArgs->numObs];

   NEW_PRINT("double", pArgs->numObs);
   pArgs->expB = new double[pArgs->numObs];

   NEW_PRINT("double", pArgs->numObs);
   pArgs->expD = new double[pArgs->numObs];

   NEW_PRINT("double", pArgs->numObs);
   pArgs->conc = new double[pArgs->numObs];

   NEW_PRINT("double", pArgs->numObs);
   pArgs->wconc = new double[pArgs->numObs];

   NEW_PRINT("double", pArgs->numObs);
   pArgs->sorb = new double[pArgs->numObs];

   NEW_PRINT("double", pArgs->numObs);
   pArgs->wsorb = new double[pArgs->numObs];
   MEM_CHECK(pArgs->wsorb);

   pTmp = strstr(pStr,"BeginLabData");
   pTmp = ISO_GetLine(pTmp, &pLine);
   xsd = sd = 1.00e00;
   for(i = 0; i < pArgs->numObs; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);

      pArgs->expA[i] = 0.00;
      pArgs->expB[i] = 1.00;
      pArgs->expD[i] = 1.00;
      if(strcmp(pType, "IndividualStdDevs") == 0)
      {
         if((pArgs->method == OREAR_METHOD) || 
            (pArgs->method == MCCAMMON_METHOD) ||
            (pArgs->method == TOTAL_ERROR_METHOD))
         {
            sscanf(pLine, "%lf %lf %lf %lf", &(pArgs->conc[i]), 
                   &(pArgs->sorb[i]), &sd, &xsd);
         }
         else if((pArgs->method == KINNIBURGH_METHOD) || (pArgs->method == ADV_KINNIBURGH_METHOD))
         {
            sscanf(pLine, "%lf %lf %lf %lf", &(pArgs->conc[i]), &sd, &(pArgs->expA[i]), &(pArgs->expB[i]));
            pArgs->expD[i] = 1.00; //not used
            pArgs->sorb[i] = 0.00; //not used
         }
         else
         {
            sscanf(pLine, "%lf %lf %lf %lf %lf %lf", &(pArgs->conc[i]), 
                   &(pArgs->sorb[i]), &sd, &(pArgs->expA[i]), 
                   &(pArgs->expB[i]), &(pArgs->expD[i]));
         }
      }
      else
      {
         sscanf(pLine, "%lf %lf %lf %lf %lf", &(pArgs->conc[i]), 
                &(pArgs->sorb[i]), &(pArgs->expA[i]), &(pArgs->expB[i]),
                &(pArgs->expD[i]));
      }

      if((pArgs->method == KINNIBURGH_METHOD) || (pArgs->method == ADV_KINNIBURGH_METHOD))
      {      
         pArgs->wconc[i] = 1.00e00/sd;
         pArgs->wsorb[i] = 1.00e00/xsd;
      }
      else
      {
         pArgs->wsorb[i] = 1.00e00/sd;
         pArgs->wconc[i] = 1.00e00/xsd;
      }
   }

   /*-------------------------------------------------------------------------------------
   Assign uniform or relative weights
   -------------------------------------------------------------------------------------*/
   if(strcmp(pType, "Uniform") == 0)
   {
      for(i = 0; i < pArgs->numObs; i++){ pArgs->wconc[i] = val;}
      for(i = 0; i < pArgs->numObs; i++){ pArgs->wsorb[i] = val;}
   }
   else if (strcmp(pType, "SorbedRelative") == 0)
   { 
      for(i = 0; i < pArgs->numObs; i++){ pArgs->wconc[i] = 1.96 / (val * pArgs->sorb[i]);}
      for(i = 0; i < pArgs->numObs; i++){ pArgs->wsorb[i] = 1.96 / (val * pArgs->sorb[i]);}
   }
   else if (strcmp(pType, "AqueousRelative") == 0)
   { 
      for(i = 0; i < pArgs->numObs; i++){ pArgs->wconc[i] = 1.96 / (convFactor * val * pArgs->conc[i]);}
      for(i = 0; i < pArgs->numObs; i++){ pArgs->wsorb[i] = 1.96 / (convFactor * val * pArgs->conc[i]);}
   }
   else if(strcmp(pType, "IndividualStdDevs") == 0)
   {
      //nothing to do, already assigned
   }
   else
   {
      LogError(ERR_BAD_ARGS, "ISO_ReadIsoFitFile() : Unknown weighting scheme, valid syntax is:");
      LogError(ERR_CONTINUE, "*****************************************************");
      LogError(ERR_CONTINUE, "Uniform <weight>");
      LogError(ERR_CONTINUE, "SorbedRelative  <relative error>");
      LogError(ERR_CONTINUE, "AqueousRelative <relative error>  <conversion factor>");
      LogError(ERR_CONTINUE, "IndividualStdDevs");
      LogError(ERR_CONTINUE, "*****************************************************");
      delete [] pStr;
      ExitProgram(1);
   }
   
   if(strstr(pStr, "PreserveOutputFiles") != NULL){pArgs->debug = 1;}
   else{pArgs->debug = 0;}

   pArgs->bHoldObs = true;
   pArgs->bHoldParams = true;
   if(strstr(pStr, "ExcludeInsensitiveParameters")   != NULL){ pArgs->bHoldParams = true;}
   if(strstr(pStr, "IncludeInsensitiveParameters")   != NULL){ pArgs->bHoldParams = false;}
   if(strstr(pStr, "ExcludeInsensitiveObservations") != NULL){ pArgs->bHoldObs = true;}
   if(strstr(pStr, "IncludeInsensitiveObservations") != NULL){ pArgs->bHoldObs = false;}

   //If debugging, write input string to file
   if(pArgs->debug == 1)
   {
      pFile = fopen("IsoFitIn.str", "w");
      i = 0;
      while(pStr[i] != (char)NULL)
      {
         fputc(pStr[i], pFile);
         i++;
      }
      fclose(pFile);
   }/* end if() */

   if(strcmp(pArgs->isoStr,  "AllIsotherms") == 0) pArgs->bFitAll = true;
   if(strcmp(pArgs->lumpStr, "yes") == 0) pArgs->bLumpedQ0 = true;

   //read in swarm parameters
   ISO_GetSwarmParams(pStr, pArgs);

   //read in solution settings, if needed
   ISO_GetSolutionSettings(pStr, pArgs);

   delete [] pStr;   
}/* end ISO_ReadIsoFitFile() */

/******************************************************************************
ISO_GetRangesSection()

Read in the Ranges section.

Returns the section contents in string format.
******************************************************************************/
char * ISO_GetRangesSection(void)
{
   char * pTmp;
   char * pLine;
   int size;
   char * pStr;
   char * pRanges;

   size = ISO_GetFileSize("IsoFitIn.txt");
   if(size <= 0)
   {
      LogError(ERR_FILE_IO, "ISO_GetRangesSection() : empty or non-existant input file");
      ExitProgram(1);
   }
   NEW_PRINT("char", size+1);
   pStr = new char[size+1];
   MEM_CHECK(pStr);
   ISO_FileToStr("IsoFitIn.txt", pStr, size);

   //check tokens
   if((pTmp=strstr(pStr, gBeginRanges)) == NULL){ delete [] pStr; return NULL;}
   if(strstr(pTmp, gEndRanges)   == NULL)
   { 
      delete [] pStr; 
      LogError(ERR_FILE_IO, "Missing token (End*Ranges)\n"); 
      printf("Missing token (End*Ranges)\n"); 
      ExitProgram(1);
   }

   pTmp = strstr(pStr, gBeginRanges);

   size = (int)strlen(pTmp);
   NEW_PRINT("char", size+1);
   pRanges = new char[size+1];
   MEM_CHECK(pRanges);
   pRanges[0] = (char)NULL;

   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strncmp(pLine, gEndRanges, strlen(gEndRanges)) != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strncmp(pLine, gEndRanges, strlen(gEndRanges)) != 0)
      {
         strcat(pRanges, pLine);
         strcat(pRanges, "\n");
      }
   }
   delete [] pStr;
   if(pRanges[0] == (char)NULL) { delete [] pRanges; return NULL;}

   return pRanges;
}/* end ISO_GetRangesSection() */

/******************************************************************************
ISO_GetFileSize()

Count the number of bytes in the file. Also sizes the maximum number of bytes
in a given line. Passing in a NULL argument will delete the data line array.
******************************************************************************/
int ISO_GetFileSize(const char * file)
{
   int maxLineSize = 0;
   int lineSize = 0;
   int fileSize = 0;
   char c;

   if(file == NULL)
   {
      delete [] gIsoLine;
      gIsoLineSize = 0;
      return 0;
   }

   FILE * pFile = fopen(file, "r");
   
   if(pFile == NULL)
   {
      printf("Couldn't open file |%s|\n", file);      
      return 0;
   }/* end if() */

   while(feof(pFile) == 0) 
   {
      fileSize++;
      c = (char)(fgetc(pFile));
      lineSize++;
      if(c == '\n')
      {
         if(lineSize > maxLineSize){
            maxLineSize = lineSize;
         }
         lineSize = 0;
      }
   }/* end while() */
   fileSize--;   

   fclose(pFile);

   //last line might not end in a carriage return
   if(lineSize > maxLineSize){
      maxLineSize = lineSize;
   }
   
   maxLineSize *= 2;
   if(maxLineSize > gIsoLineSize){
      delete [] gIsoLine;
      gIsoLineSize = maxLineSize;
      NEW_PRINT("char *", gIsoLineSize);
      gIsoLine = new char[gIsoLineSize];
      MEM_CHECK(gIsoLine);
   }

   return (fileSize);
}/* end ISO_GetFileSize() */

/******************************************************************************
ISO_FileToStr()

read file into pStr, pStr must be pre-allocated
******************************************************************************/
void ISO_FileToStr(const char * file, char * pStr, int size)
{
   bool skip_line = false;
   int i, j;
   char c;
   FILE * pFile = fopen(file, "r");
   if(pFile == NULL)
   {
      printf("Couldn't open file |%s|\n", file);
      pStr[0] = (char)NULL;
      return;
   }/* end if() */

   //copy line from file to string, omit comments
   j = 0;
   for(i = 0; i < size; i++)
   {
      c = fgetc(pFile);
      if(c == '#'){ skip_line = true;}
      if((c == '\n') || (c == '\r')){ skip_line = false;}
      if(skip_line == false)
      {      
         pStr[j] = c;
         j++;
      }
   }
   pStr[j] = (char)NULL;
   fclose(pFile);   
}/* end FileToStr() */

/******************************************************************************
ISO_GetLine()

read next line of string and store in the global array (gIsoLine), which is
returned via the pLine pointer reference. Return value is a pointer to the
string after it has been advanced past the line.
******************************************************************************/
char * ISO_GetLine(char * pStr, char ** pLine)
{
   char * pTmp = gIsoLine;
   *pTmp = (char)NULL;

   if(pLine != NULL){
      *pLine = pTmp;
   }
   
   while((*pStr != (char)NULL) && (*pStr != '\n') && (*pStr != '\r'))
   {
      if(pTmp != NULL){ *pTmp = *pStr; pTmp++;}
      pStr++;
   }
   if(pTmp != NULL) {*pTmp = (char)NULL;}

   if(*pStr == (char)NULL){ return pStr;}
   while((*pStr == '\n') || (*pStr == '\r'))
   {
      pStr++;
      if(*pStr == (char)NULL){ return pStr;}
   }
   return pStr;
}/* end GetLine() */

/******************************************************************************
ISO_CreateTemplateFile()

create a template file
******************************************************************************/
void ISO_CreateTemplateFile(IsoGlobStruct * pArgs)
{
   int i;
   FILE * pFile = fopen(ISO_TPL_FILE, "w");
   if(pFile == NULL){ FileOpenFailure("ISO_CreateTemplateFile()", ISO_TPL_FILE);}

   if(pArgs->bLumpedQ0 == true){
      fprintf(pFile, "Isotherm Template File, AutoGenerated by Ostrich\n" \
                  "IsothermType %s\n" \
                  "KinniburghLossTerm  XVal\n\n" \
                  "BeginLinearIsotherm\nKd KdVal\nEndLinearIsotherm\n\n" \
                  "BeginLangmuirIsotherm\nb*Q0 b*Q0Val\nb  bVal\nEndLangmuirIsotherm\n\n" \
                  "BeginFreundlichIsotherm\nKf KfVal\n(1/n)  (1/n)Val\nEndFreundlichIsotherm\n\n" \
                  "BeginPolanyi-PartitionIsotherm\nQ0 Q0Val\na  aVal\nb  bVal\nKp KpVal\nSw %s\nEndPolanyi-PartitionIsotherm\n\n" \
                  "BeginLangmuir-PartitionIsotherm\nb*Q0 b*Q0Val\nb  bVal\nKp KpVal\nEndLangmuir-PartitionIsotherm\n\n" \
                  "BeginPolanyiIsotherm\nQ0 Q0Val\na  aVal\nb  bVal\nSw %s\nEndPolanyiIsotherm\n\n" \
                  "BeginBET_Isotherm\nb*Q0 b*Q0Val\nb  bVal\nSw %s\nEndBET_Isotherm\n\n" \
                  "BeginTothIsotherm\nb*Q0 b*Q0Val\nb  bVal\nn  nVal\nEndTothIsotherm\n\n" \
                  "BeginLangmuir-FreundlichIsotherm\nQ0 Q0Val\nb  bVal\n(1/n)  (1/n)Val\nEndLangmuir-FreundlichIsotherm\n\n" \
                  "BeginFreundlich-PartitionIsotherm\nKf KfVal\nKp  KpVal\n(1/n)  (1/n)Val\nEndFreundlich-PartitionIsotherm\n\n" \
                  "BeginDualLangmuirIsotherm\nb1*Q01 b1*Q01Val\nb1  b1Val\nb2*Q02 b2*Q02Val\nb2  b2Val\nEndDualLangmuirIsotherm\n\n" \
                  "BeginOrearIsotherm\na aVal\nb  bVal\nEndOrearIsotherm\n\n" \
                  "BeginMcCammonIsotherm\nA AVal\nB  BVal\nC   CVal\n_E_   EVal\nF   FVal\nEndMcCammonIsotherm\n\n" \
                  "Name	C\nBeginConcentrations\n", pArgs->isoStr, pArgs->solStr, pArgs->solStr, pArgs->solStr);}
   else{
      fprintf(pFile, "Isotherm Template File, AutoGenerated by Ostrich\n" \
                  "IsothermType %s\n" \
                  "KinniburghLossTerm  XVal\n\n" \
                  "BeginLinearIsotherm\nKd KdVal\nEndLinearIsotherm\n\n" \
                  "BeginLangmuirIsotherm\nQ0 Q0Val\nb  bVal\nEndLangmuirIsotherm\n\n" \
                  "BeginFreundlichIsotherm\nKf KfVal\n(1/n)  (1/n)Val\nEndFreundlichIsotherm\n\n" \
                  "BeginPolanyi-PartitionIsotherm\nQ0 Q0Val\na  aVal\nb  bVal\nKp KpVal\nSw %s\nEndPolanyi-PartitionIsotherm\n\n" \
                  "BeginLangmuir-PartitionIsotherm\nQ0 Q0Val\nb  bVal\nKp KpVal\nEndLangmuir-PartitionIsotherm\n\n" \
                  "BeginPolanyiIsotherm\nQ0 Q0Val\na  aVal\nb  bVal\nSw %s\nEndPolanyiIsotherm\n\n" \
                  "BeginBET_Isotherm\nQ0 Q0Val\nb  bVal\nSw %s\nEndBET_Isotherm\n\n" \
                  "BeginTothIsotherm\nQ0 Q0Val\nb  bVal\nn  nVal\nEndTothIsotherm\n\n" \
                  "BeginLangmuir-FreundlichIsotherm\nQ0 Q0Val\nb  bVal\n(1/n)  (1/n)Val\nEndLangmuir-FreundlichIsotherm\n\n" \
                  "BeginFreundlich-PartitionIsotherm\nKf KfVal\nKp  KpVal\n(1/n)  (1/n)Val\nEndFreundlich-PartitionIsotherm\n\n" \
                  "BeginDualLangmuirIsotherm\nQ01 Q01Val\nb1  b1Val\nQ02 Q02Val\nb2  b2Val\nEndDualLangmuirIsotherm\n\n" \
                  "BeginOrearIsotherm\na aVal\nb  bVal\nEndOrearIsotherm\n\n" \
                  "BeginMcCammonIsotherm\nA AVal\nB  BVal\nC   CVal\n_E_   EVal\nF   FVal\nEndMcCammonIsotherm\n\n" \
                  "Name	C\nBeginConcentrations\n", pArgs->isoStr, pArgs->solStr, pArgs->solStr, pArgs->solStr);}

   if((pArgs->method == ISOTHERM_METHOD) || (pArgs->method == KINNIBURGH_METHOD) || 
                                            (pArgs->method == ADV_KINNIBURGH_METHOD))
   {
      for(i = 0; i < pArgs->numObs; i++){ fprintf(pFile,"obs%d\t%E\n", i, pArgs->conc[i]);}
   }
   else if(pArgs->method == TOTAL_ERROR_METHOD) //treat concentrations as Ostrich parameters
   {
      for(i = 0; i < pArgs->numObs; i++){ fprintf(pFile,"obs%d\tConc%d_Val\n", i, i);}
   }
   else //Orear/McCammon requires additional concentration information
   {
      for(i = 0; i < pArgs->numObs; i++)
      { 
         fprintf(pFile,"obs%d\t%E\t%E\t%E\t%E\n", i, pArgs->conc[i], 
                 pArgs->sorb[i], pArgs->wconc[i], pArgs->wsorb[i]);
      }
   }

   fprintf(pFile, "EndConcentrations\n");

   if(pArgs->method == KINNIBURGH_METHOD)
   {
      fprintf(pFile, "\nBeginKinniburghMethod\n");
      fprintf(pFile, "MaxBisections  %d\n", pArgs->maxBisections);
      fprintf(pFile, "EndKinniburghMethod\n");

      fprintf(pFile, "\nBeginExperimentalConstants\n");
      for(i = 0; i < pArgs->numObs; i++){ fprintf(pFile,"%E  %E  %E\n", pArgs->expA[i], pArgs->expB[i], pArgs->expD[i]);}
      fprintf(pFile, "EndExperimentalConstants\n");
   }
   else if(pArgs->method == ADV_KINNIBURGH_METHOD)
   {
      fprintf(pFile, "\nBeginAdvancedKinniburghMethod\n");
      fprintf(pFile, "MaxBisections  %d\n", pArgs->maxBisections);
      fprintf(pFile, "EndAdvancedKinniburghMethod\n");

      fprintf(pFile, "\nBeginExperimentalConstants\n");
      for(i = 0; i < pArgs->numObs; i++){ fprintf(pFile,"%E  %E  %E\n", pArgs->expA[i], pArgs->expB[i], pArgs->expD[i]);}
      fprintf(pFile, "EndExperimentalConstants\n");
   }
   else if(pArgs->method == OREAR_METHOD)
   {
      fprintf(pFile, "\nBeginOrearMethod\n");
      fprintf(pFile, "MaxBisections  %d\n", pArgs->maxBisections);
      fprintf(pFile, "EndOrearMethod\n");
   }
   else if(pArgs->method == MCCAMMON_METHOD)
   {
      fprintf(pFile, "\nBeginMcCammonMethod\n");
      fprintf(pFile, "MaxBisections  %d\n", pArgs->maxBisections);
      fprintf(pFile, "EndMcCammonMethod\n");
   }

   fclose(pFile);
}/* end ISO_CreateTemplateFile() */

/******************************************************************************
ISO_CreateOstrichFile()

create a Ostrich input file
******************************************************************************/
void ISO_CreateOstrichFile(const char * pProgType, IsoParamList * pList, bool stats, 
                           IsoGlobStruct * pArgs)
{
   int i, j, col, xcol, np = 0;
   IsoParamList * pTmp;
   double * pObs, * pExtra;
   double * pWeight, * pXwght;

   FILE * pFile = fopen(ISO_OSTIN_FILE, "w");

   for(pTmp = pList; pTmp != NULL; pTmp = pTmp->pNxt){ np++;}
   //set PSO defaults unless user has already set up
   if(pArgs->popSize == 0){ pArgs->popSize = 20*np;}
   if(pArgs->maxGens == 0){ pArgs->maxGens = 20*np;}

   if(pFile == NULL){ FileOpenFailure("ISO_CreateOstrichFile()", ISO_OSTIN_FILE);}
   fprintf(pFile, "#Configuration File for Ostrich Program\n\n" \
                  "ProgramType %s\n\n" \
                  "ModelSubdir    .\n\n" \
                  "NumDigitsOfPrecision 16\n\n" \
                  "BeginFilePairs\n%s\t%s\nEndFilePairs\n\n", pProgType, ISO_TPL_FILE, ISO_IN_FILE);

   if(stats == true)
   {
      fprintf(pFile, "CheckSensitivities no\n");
   }

   if(pArgs->method == ISOTHERM_METHOD)
   {
      fprintf(pFile, "ModelExecutable    Isotherm()\n\n");
      pObs = pArgs->sorb;
      pWeight = pArgs->wsorb;
      pExtra = NULL;
      pXwght = NULL;
      col = 3;
   }
   else if(pArgs->method == TOTAL_ERROR_METHOD)
   {
      fprintf(pFile, "ModelExecutable    Isotherm()\n\n");
      pObs = pArgs->sorb;
      pWeight = pArgs->wsorb;
      pExtra = pArgs->conc;
      pXwght = pArgs->wconc;
      col = 3;
      xcol = 2;
   }
   else if (pArgs->method == OREAR_METHOD)
   {
      fprintf(pFile, "ModelExecutable    Orear()\n\n");
      pObs = pArgs->sorb;
      pWeight = pArgs->wsorb;
      pExtra = pArgs->conc;
      pXwght = pArgs->wconc;
      col = 3;
      xcol = 2;
   }
   else if (pArgs->method == MCCAMMON_METHOD)
   {
      fprintf(pFile, "ModelExecutable    McCammon()\n\n");
      pObs = pArgs->sorb;
      pWeight = pArgs->wsorb;
      pExtra = pArgs->conc;
      pXwght = pArgs->wconc;
      col = 3;
      xcol = 2;
   }
   else if (pArgs->method == KINNIBURGH_METHOD)
   {
      fprintf(pFile, "ModelExecutable    Kinniburgh()\n\n");
      pObs = pArgs->conc;
      pWeight = pArgs->wconc;
      pExtra = NULL;
      pXwght = NULL;
      col = 2;
   }
   else if (pArgs->method == ADV_KINNIBURGH_METHOD)
   {
      fprintf(pFile, "ModelExecutable    AdvancedKinniburgh()\n\n");
      pObs = pArgs->conc;
      pWeight = pArgs->wconc;
      pExtra = NULL;
      pXwght = NULL;
      col = 2;
   }
   else //default to standard method
   {
      fprintf(pFile, "ModelExecutable    Isotherm()\n\n");
      pObs = pArgs->sorb;
      pWeight = pArgs->wsorb;
      pExtra = NULL;
      pXwght = NULL;
      col = 3;
   }
   fprintf(pFile, "#Parameter Specification\n" \
                  "BeginParams\n" \
                  "#parameter	init.	low	high	tx_in  tx_ost	tx_out\n");

   pTmp = pList;
   while(pTmp != NULL)
   {
      fprintf(pFile,"%s\t%E\t%E\t%E\t%s\t%s\t%s\n", pTmp->name, pTmp->init, 
              pTmp->lwr, pTmp->upr, pTmp->txin, pTmp->txost, pTmp->txout);
      pTmp = pTmp->pNxt;
   }/* end while() */

   fprintf(pFile,"EndParams\n\n");

   /* -------------------------------------------------------------------
   TotalError Method treats aqueous concentrations as Ostrich parameters
   Seed the PSO with observed C as an initial quess.
   --------------------------------------------------------------------- */
   if(pArgs->method == TOTAL_ERROR_METHOD)
   {
      fprintf(pFile,"BeginInitParams\n");
      pTmp = pList;
      while(pTmp != NULL)
      {
         fprintf(pFile,"%E  ", pTmp->init);
         pTmp = pTmp->pNxt;
      }/* end while() */
      fprintf(pFile,"\nEndInitParams\n\n");
   }/* end if() */

   fprintf(pFile,"#Observation Configuration\n" \
                 "BeginObservations\n" \
                  "#observation	value	weight	file		keyword		line	column\n");

   for(i = 0; i < pArgs->numObs; i++)
   {
      fprintf(pFile, "obs%d\t%E\t%E\t%s\tConcentration\t%d\t%d\n", 
                     i, pObs[i], pWeight[i], ISO_OUT_FILE, (i+1), col);
   }/* end for() */

   //Orear, McCammon and TotalError consider WSSE of both q and C
   if(pExtra != NULL)
   {
      for(i = 0; i < pArgs->numObs; i++)
      {
         j = i + pArgs->numObs;
         fprintf(pFile, "obs%d\t%E\t%E\t%s\tConcentration\t%d\t%d\n", 
                        j, pExtra[i], pXwght[i], ISO_OUT_FILE, (i+1), xcol);
      }/* end for() */
   }/* end if() */

   fprintf(pFile,
"EndObservations\n\n\
#Configuration for Levenberg-Marquardt algorithm\n\
BeginLevMar\n\
  InitialLambda    10.0\n\
  LambdaScaleFactor    1.1\n\
  MoveLimit    0.1\n\
  AlgorithmConvergenceValue    1E-10\n\
  LambdaPhiRatio    0.3\n\
  LambdaRelReduction    0.01\n\
  MaxLambdas    10\n\
  MaxIterations    100\n\
EndLevMar\n\n\
\
BeginParticleSwarm\n\
  SwarmSize %d\n\
  InertiaReductionRate linear\n\
  NumGenerations %d\n\
  InitPopulationMethod LHS\n\
  ConvergenceVal -1.00\n\
  EndParticleSwarm\n\n", pArgs->popSize, pArgs->maxGens);

   fprintf(pFile,"#Powell's Algorithm Configuration\n" \
                 "BeginPowellAlg\n" \
                 "ConvergenceVal 1E-10\n" \
                 "MaxIterations 200\n" \
                 "EndPowellAlg\n\n" \
                 "#Configuration of One-Dimensional Search\n" \
                 "Begin1dSearch\n" \
                 "1dSearchConvergeVal 1.000000E-006\n" \
                 "1dSearchMethod Brent\n" \
                 "End1dSearch\n\n");

   fprintf(pFile,"BeginMathAndStats\nDiffType    forward\nDiffIncType    value-relative\nDiffIncrement    ");
   //fprintf(pFile,"BeginMathAndStats\nDiffType    forward\nDiffRelIncrement    ");
   pTmp = pList;
   while(pTmp != NULL)
   {
      fprintf(pFile,"%E  ", pTmp->step);
      pTmp = pTmp->pNxt;
   }/* end while() */
   fprintf(pFile,"\n");

   if(stats == true)
   {
      fprintf(pFile,"StdDev\nStdErr\nCorrCoeff\nBeale\nLinssen\nCooksD\nDFBETAS\nConfidence\nNormPlot\n" \
                    "Sensitivity\nMatrices\nRunsTest\nAutorunFunction\nBestBoxCox\nMMRI\n");
      if(pArgs->bHoldObs == true)
      {
         fprintf(pFile,"ExcludeInsensitiveObservations\n");
      }
      else
      {
         fprintf(pFile,"IncludeInsensitiveObservations\n");
      }
      if(pArgs->bHoldParams == true)
      {
         fprintf(pFile,"ExcludeInsensitiveParameters\n");
      }
      else
      {
         fprintf(pFile,"IncludeInsensitiveParameters\n");
      }
   }
   fprintf(pFile,"EndMathAndStats\n");

   // append any user-defined extras (e.g. random seed)
   char buf[1000];
   FILE * pExtras = fopen("OstInExtras.txt", "r");
   if(pExtras != NULL)
   {
      while(!feof(pExtras))
      {
        fgets(buf, 1000, pExtras);
        fputs(buf, pFile);
        buf[0]=NULLSTR;
      }
      fclose(pExtras);
   }/* end if() */
   fclose(pFile);
}/* end ISO_CreateOstrichFile() */

/******************************************************************************
ISO_CreateParamList()

fill in the parameter list with values appropriate for the given Isotherm
******************************************************************************/
void ISO_CreateParamList(IsoParamList * pList, IsoGlobStruct * pArgs)
{
   IsoParamList * pTmp;
   char * pCur;

   if(strcmp(pArgs->isoStr, "LinearIsotherm") == 0)
   {
      strcpy(pList->name, "KdVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;
      pList->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "LangmuirIsotherm") == 0)
   {
      if(pArgs->bLumpedQ0 == true){ strcpy(pList->name, "b*Q0Val");}
      else                 { strcpy(pList->name, "Q0Val");}
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 50.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "FreundlichIsotherm") == 0)
   {
      printf("Adding Kf val\n");

      strcpy(pList->name, "KfVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      printf("Adding (1/n) val\n");

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "(1/n)Val");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.50;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+0;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "DualLangmuirIsotherm") == 0)
   {
      if(pArgs->bLumpedQ0 == true){ strcpy(pList->name, "b1*Q01Val");}
      else                 { strcpy(pList->name, "Q01Val");}
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "b1Val");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 50.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      if(pArgs->bLumpedQ0 == true){ strcpy(pTmp->name, "b2*Q02Val");}
      else                 { strcpy(pTmp->name, "Q02Val");}
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 100.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "b2Val");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 50.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "Polanyi-PartitionIsotherm") == 0)
   {
      strcpy(pList->name, "KpVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "Q0Val");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 100.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "aVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.10;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+0;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 2.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+1;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "Langmuir-PartitionIsotherm") == 0)
   {
      strcpy(pList->name, "KpVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      if(pArgs->bLumpedQ0 == true){ strcpy(pTmp->name, "b*Q0Val");}
      else                 { strcpy(pTmp->name, "Q0Val");}

      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 50.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 100.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "BET_Isotherm") == 0)
   {
      if(pArgs->bLumpedQ0 == true){ strcpy(pList->name, "b*Q0Val");}
      else                 { strcpy(pList->name, "Q0Val");}
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;

      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 100.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e+0;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "TothIsotherm") == 0)
   {
      if(pArgs->bLumpedQ0 == true){ strcpy(pList->name, "b*Q0Val");}
      else                 { strcpy(pList->name, "Q0Val");}
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;
      pTmp = pList->pNxt;

      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 50.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "nVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.50;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+0;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "Langmuir-FreundlichIsotherm") == 0)
   {
      strcpy(pList->name, "Q0Val");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;
      pTmp = pList->pNxt;

      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 50.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "(1/n)Val");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.50;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+0;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "PolanyiIsotherm") == 0)
   {
      strcpy(pList->name, "Q0Val");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;
      pTmp = pList->pNxt;

      strcpy(pTmp->name, "aVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.10;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+0;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;
      pTmp = pTmp->pNxt;

      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 2.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+1;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "Freundlich-PartitionIsotherm") == 0)
   {
      strcpy(pList->name, "KfVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "log10");
      pList->init = 100.00;
      pList->step = DEFAULT_STEP;
      pList->upr = 1e+6;
      pList->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "(1/n)Val");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.50;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+0;
      pTmp->lwr = 1e-6;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;

      pTmp = pTmp->pNxt;
      strcpy(pTmp->name, "KpVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 100.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1e+6;
      pTmp->lwr = 1e-6;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "OrearIsotherm") == 0)
   {
      strcpy(pList->name, "aVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "none");
      pList->init = 0.50;
      pList->step = DEFAULT_STEP;
      pList->upr = 1.00;
      pList->lwr = 0.00;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "bVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "log10");
      pTmp->init = 100.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 1.00e+6;
      pTmp->lwr = 1.00e+0;
      pTmp->pNxt = NULL;
   }
   else if(strcmp(pArgs->isoStr, "McCammonIsotherm") == 0)
   {
      strcpy(pList->name, "AVal");
      strcpy(pList->txin, "none");
      strcpy(pList->txout, "none");
      strcpy(pList->txost, "none");
      pList->init = 0.0125;
      pList->step = DEFAULT_STEP;
      pList->upr = 0.025;
      pList->lwr = 0.00;

      NEW_PRINT("IsoParamList", 1);
      pList->pNxt = new IsoParamList;

      pTmp = pList->pNxt;
      strcpy(pTmp->name, "BVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = -0.25;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 0.00;
      pTmp->lwr = -0.50;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;

      pTmp = pTmp->pNxt;
      strcpy(pTmp->name, "CVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.0125;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 0.025;
      pTmp->lwr = 0.00;

      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;

      pTmp = pTmp->pNxt;
      strcpy(pTmp->name, "EVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = -0.25;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 0.00;
      pTmp->lwr = -0.50;
      pTmp->pNxt = NULL;
   }
   else
   {
      LogError(ERR_BAD_ARGS, "ISO_CreateParamList() : Unknown Isotherm type");
      printf("ISO_CreateParamList() : Unknown Isotherm type");
      ISO_DestroyIsoParamList(pList);
      ExitProgram(1);
   }

   //advanced Kinniburgh method has additional loss term parameter
   if(pArgs->method == ADV_KINNIBURGH_METHOD)
   {
      NEW_PRINT("IsoParamList", 1);
      pTmp->pNxt = new IsoParamList;

      pTmp = pTmp->pNxt;
      strcpy(pTmp->name, "XVal");
      strcpy(pTmp->txin, "none");
      strcpy(pTmp->txout, "none");
      strcpy(pTmp->txost, "none");
      pTmp->init = 0.00;
      pTmp->step = DEFAULT_STEP;
      pTmp->upr = 0.00;
      pTmp->lwr = 1.00;
      pTmp->pNxt = NULL;
   }

   /*------------------------------------------------------
   Open up input file and search through it for the 
   inclusion of a Ranges section, which will have lower 
   and upper bounds for parameters and also transformations,
   if requested.
   -------------------------------------------------------*/
   char * pRanges = ISO_GetRangesSection();
   if(pRanges != NULL)
   {
      for(pCur = pRanges; *pCur != (char)NULL; pCur = ISO_GetLine(pCur, NULL))
      {
         //step through param list, matching names with ranges section
         for(pTmp = pList; pTmp != NULL; pTmp = pTmp->pNxt)
         {
            int len = (int)strlen(pTmp->name)-3;
            if(strncmp(pCur, pTmp->name, len) == 0)
            {
               pCur += len;
               sscanf(pCur, "%lf %lf %s %lf", &(pTmp->lwr), &(pTmp->upr), pTmp->txost, &(pTmp->step));
               pTmp->init = 0.5*(pTmp->lwr + pTmp->upr);
               //if(strstr(pCur, "log") != NULL){ strcpy(pTmp->txost, "log10");}
               //else                           { strcpy(pTmp->txost, "none");}
            }
         }
      }
      delete [] pRanges;
   }

   /* ---------------------------------------------------------------------
   TotalError Method treats aqueous concentrations as Ostrich parameters
   --------------------------------------------------------------------- */
   if(pArgs->method == TOTAL_ERROR_METHOD)
   {
      int i;
      double init;

      //locate end of list
      for(pTmp = pList; pTmp->pNxt != NULL; pTmp = pTmp->pNxt);

      for(i = 0; i < pArgs->numObs; i++)
      {
         init = pArgs->conc[i] + 1E-10;

         NEW_PRINT("IsoParamList", 1);
         pTmp->pNxt = new IsoParamList;

         pTmp = pTmp->pNxt;
         sprintf(pTmp->name, "Conc%d_Val", i);
         strcpy(pTmp->txin, "none");
         strcpy(pTmp->txout, "none");
         strcpy(pTmp->txost, "none");
         pTmp->init = init;
         pTmp->step = DEFAULT_STEP;
         pTmp->upr = 2.0*init;
         pTmp->lwr = 0.5*init;
         pTmp->pNxt = NULL;         
      }/* end for() */
   }/* end if() */
}/* end ISO_CreateParamList() */

/******************************************************************************
ISO_RefreshParamList()

Read optimal parameters from Ostrich output file. Also switch from log10 to 
none on the Ostrich transformation.
******************************************************************************/
void ISO_RefreshParamList(IsoParamList * pList)
{
   int j;
   IsoParamList * pCur;
   char * pTmp, * pTok;
   char varStr[DEF_STR_SZ];
   char * pLine;
   char valStr[DEF_STR_SZ];
   int size;
   char * pStr;

   size = ISO_GetFileSize(ISO_OSTOUT_FILE);
   if(size <= 0)
   {
      LogError(ERR_FILE_IO, "ISO_RefreshParamList() : empty or nonexistent output file");
      ExitProgram(1);
   }
   NEW_PRINT("char", size+1);
   pStr = new char[size+1];
   MEM_CHECK(pStr);
   ISO_FileToStr(ISO_OSTOUT_FILE, pStr, size);

   //locate "Optimal Parameter Set"
   pTmp = strstr(pStr, "Optimal Parameter Set");
   if(pTmp == NULL)
   {
      LogError(ERR_FILE_IO, "ISO_RefreshParamList() : couldn't locate Optimal Parameter Set");
      delete [] pStr;
      ExitProgram(1);
   }

   //read in each optimal parameter
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strstr(pLine, "Observation Residuals") == NULL)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(*pTmp == (char)NULL)
      {
         LogError(ERR_FILE_IO, "ISO_RefreshParamList() : couldn't locate Observation Residuals");
         delete [] pStr;
         ExitProgram(1);
      }/* end if() */
      if(strstr(pLine, ":") != NULL)
      {
         //extract variable name and value
         pTok = pLine;
         j = ExtractColString(pTok, varStr, ':');
         j = ValidateExtraction(j, 1, 1, "ISO_RefreshParamList()");
         pTok += j;
         ExtractColString(pTok, valStr, ':');
         MyTrim(varStr);
         MyTrim(valStr);

         pCur = pList;
         while(pCur != NULL)
         {
            if(strcmp(pCur->name, varStr) == 0)
            {
               pCur->init = atof(valStr);
               strcpy(pCur->txost, "none");
               break;
            }
            pCur = pCur->pNxt;
         }/* end while() */
      }/* end if() */
   }/* end while() */
}/* end ISO_RefreshParamList() */

/******************************************************************************
ISO_DestroyIsoParamList()

Recursively free up the parameter list.
******************************************************************************/
void ISO_DestroyIsoParamList(IsoParamList * pList)
{
   if(pList == NULL){ return;}
   if(pList->pNxt != NULL){ ISO_DestroyIsoParamList(pList->pNxt);}
   delete pList;
}/* end ISO_DestroyIsoParamList() */

/******************************************************************************
ISO_GetSolutionSettings()

Read in solution settings. These are parameters that pertain to the selected 
solution method (Orear, McCammon, Kinniburgh, or AdvancedKinniburgh)
******************************************************************************/
void ISO_GetSolutionSettings(char * pStr, IsoGlobStruct * pArgs)
{
   char * pTmp = NULL;   
   char * pLine;
   bool bRatio = false, bVolume = false, bMax = false;

   if(pArgs->method == ISOTHERM_METHOD){ return;}
   if(pArgs->method == TOTAL_ERROR_METHOD){ return;}

   /*----------------------------------------------
   Read in common parameters
   ----------------------------------------------*/
   pTmp = pStr;
   pArgs->maxBisections = 50;
   do{ 
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strncmp(pLine, "MaxBisections", 13) == 0)
      {
         pArgs->maxBisections = atoi(&(pLine[13]));
         bMax = true;
      }
      if(bMax == true) break;
   }while(*pTmp != (char)NULL);
}/* end ISO_GetSolutionSettings() */
