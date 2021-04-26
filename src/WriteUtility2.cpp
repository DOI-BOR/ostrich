/******************************************************************************
File      : WriteUtility.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

This file provides a unifying interface for the various algorithms to write output
both to file and to stdout.

Version History
03-08-04    lsm   created
07-08-04    lsm   added objective function output, tied parameter output and
                  added date and time of build to version string
08-11-04    lsm   upped version, added date and time of build, added PATO and
                  algorithm metrics support
01-18-05    lsm   upped version, added support for GCOP
01-01-07    lsm   upped version, added support for ModelABC
06-22-07    lsm   added dual-Langmuir isotherm, upped version to 1.3
******************************************************************************/
#include <mpi.h>
#include <string.h>

#include "Algorithm.h"
#include "ParameterABC.h"
#include "PumpAndTreat.h"
#include "ResponseVarGroup.h"
#include "GenConstrainedOpt.h"
#include "ObservationGroup.h"
#include "ParameterGroup.h"

#include "WriteUtility.h"
#include "SuperMuseUtility.h"
#include "Utility.h"
#include "Exception.h"

//PATO output routines
void WriteCostToFile2(FILE * pFile, Algorithm *algo);
void WriteConstraintsToFile2(FILE * pFile, Algorithm *algo);
void WriteWellsToFile2(FILE * pFile, Algorithm *algo);

/******************************************************************************
WriteDisclaimer()

Write standard GPL/FSF disclaimer.
******************************************************************************/
void WriteDisclaimer2(FILE * pFile)
{
   int day = 0, month = 0, year = 0;
   char monthStr[10];

   sscanf(__DATE__, "%s %d %d", monthStr, &day, &year);
   if(strcmp(monthStr, "Jan") == 0) month = 1;
   else if(strcmp(monthStr, "Jan") == 0) month = 1;
   else if(strcmp(monthStr, "Feb") == 0) month = 2;
   else if(strcmp(monthStr, "Mar") == 0) month = 3;
   else if(strcmp(monthStr, "Apr") == 0) month = 4;
   else if(strcmp(monthStr, "May") == 0) month = 5;
   else if(strcmp(monthStr, "Jun") == 0) month = 6;
   else if(strcmp(monthStr, "Jul") == 0) month = 7;
   else if(strcmp(monthStr, "Aug") == 0) month = 8;
   else if(strcmp(monthStr, "Sep") == 0) month = 9;
   else if(strcmp(monthStr, "Oct") == 0) month = 10;
   else if(strcmp(monthStr, "Nov") == 0) month = 11;
   else if(strcmp(monthStr, "Dec") == 0) month = 12;

#ifdef ISOFIT_BUILD
fprintf(pFile,
"--------------------------------------------------------------------------\n \
ISOFIT version %02d.%02d.%02d (Built %s @ %s)\n\n \
A computer program for isotherm fitting.\n\n \
Author             L. Shawn Matott\n \
Copyright (C) 2006 L. Shawn Matott\n\n \
This program is free software; you can redistribute \n \
it and/or modify it under the terms of the GNU  \n \
General Public License as published by the Free \n \
Software Foundation; either version 2 of the \n \
License, or(at your option) any later version. \n\n \
This program is distributed in the hope that it will \n \
be useful, but WITHOUT ANY WARRANTY; without even \n \
the implied warranty of MERCHANTABILITY or FITNESS \n \
FOR A PARTICULAR PURPOSE. See the GNU General Public \n \
License for more details. \n\n \
You should have received a copy of the GNU General \n \
Public License along with this program; if not, \n \
write to the Free Software Foundation, Inc., 59 \n \
Temple Place, Suite 330, Boston, MA 02111-1307 USA \n\
--------------------------------------------------------------------------\n\n",
year-2000,month, day, __DATE__, __TIME__);
#endif

fprintf(pFile,
"--------------------------------------------------------------------------\n \
OSTRICH version %02d.%02d.%02d (Built %s @ %s)\n\n \
A computer program for model-independent calibration and optimization.\n\n \
Author             L. Shawn Matott\n \
Copyright (C) 2007 L. Shawn Matott\n\n \
This program is free software; you can redistribute \n \
it and/or modify it under the terms of the GNU  \n \
General Public License as published by the Free \n \
Software Foundation; either version 2 of the \n \
License, or(at your option) any later version. \n\n \
This program is distributed in the hope that it will \n \
be useful, but WITHOUT ANY WARRANTY; without even \n \
the implied warranty of MERCHANTABILITY or FITNESS \n \
FOR A PARTICULAR PURPOSE. See the GNU General Public \n \
License for more details. \n\n \
You should have received a copy of the GNU General \n \
Public License along with this program; if not, \n \
write to the Free Software Foundation, Inc., 59 \n \
Temple Place, Suite 330, Boston, MA 02111-1307 USA \n\
--------------------------------------------------------------------------\n\n",
year-2000,month, day, __DATE__, __TIME__);
}/* end writeDisclaimer() */

/******************************************************************************
WriteSetup()

Write out Ostrich Setup
******************************************************************************/
void WriteSetup2(Algorithm *algo, const char * algStr) {
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);
   remove(fileName);

   //write out disclaimers
   pFile = fopen(fileName, "a");
   WriteDisclaimer2(pFile);
   WriteSetupToFile2(pFile, algo, algStr);
   WriteSuperMuseSetupToFile(pFile);
   fclose(pFile);

   WriteDisclaimer2(stdout);
   WriteSetupToFile2(stdout, algo, algStr);
   WriteSuperMuseSetupToFile(stdout);
}/* end WriteSetup() */

/******************************************************************************
WriteSetupNoDisclaimer()

Write out Ostrich Setup without disclaimer.
******************************************************************************/
void WriteSetupNoDisclaimer2(Algorithm *algo, const char * algStr) {
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);

   //write out disclaimers
   pFile = fopen(fileName, "a");
   WriteSetupToFile2(stdout, algo, algStr);
   fclose(pFile);

   //WriteSetupToFile(stdout, algo, algStr);
}/* end WriteSetupNoDisclaimer() */

/******************************************************************************
WriteSetupToFile()

Write out Ostrich Setup
******************************************************************************/
void WriteSetupToFile2(FILE * pFile, Algorithm *algo, const char * algStr)
{
      fprintf(pFile, "Ostrich Setup\n");
      fprintf(pFile, "Model                  : %s\n", algo->GetModelStr());
      fprintf(pFile, "Algorithm              : %s\n", algStr);	  
      fprintf(pFile, "Objective Function     : %s\n", algo->GetObjFuncStr());
      fprintf(pFile, "Number of Parameters   : %d\n", algo->GetParamGroupPtr()->GetNumParams());
      fprintf(pFile, "Number of Tied Params  : %d\n", algo->GetParamGroupPtr()->GetNumTiedParams());
      fprintf(pFile, "Number of Observations : ");
      if(algo->GetObsGroupPtr() == NULL){fprintf(pFile, "0\n");}
      else {fprintf(pFile, "%d\n", algo->GetObsGroupPtr()->GetNumObs());}
	   fprintf(pFile, "Seed for Random Nums.  : %u\n", GetRandomSeed());
       algo->GetObjFuncPtr()->WriteSetupToFile(pFile);
      fprintf(pFile, "\n");
}/* end WriteSetupToFile() */

/******************************************************************************
WriteBanner()

Write out iteration banner
******************************************************************************/
void WriteBanner2(Algorithm *algo, const char * pBef, const char * pAft)
{
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);

   pFile = fopen(fileName, "a");
   WriteBannerToFile2(pFile, algo, pBef, pAft);
   fclose(pFile);

   WriteBannerToFile2(stdout, algo, pBef, pAft);
}/* end WriteBanner() */

/******************************************************************************
WriteBannerToFile()

Write out iteration banner
******************************************************************************/
void WriteBannerToFile2(FILE * pFile, Algorithm *algo, const char * pBef, const char * pAft)
{
   ObservationGroup * pObsGroup;
   ResponseVarGroup * pRespVarGroup;
   ObjectiveFunction * pObjFunc;

   pObjFunc = algo->GetObjFuncPtr();
   pObsGroup = algo->GetObsGroupPtr();
   pRespVarGroup = NULL;
   if(pObjFunc != NULL)
      pRespVarGroup = (ResponseVarGroup *)(pObjFunc->GetResponseVarGroup());


   fprintf(pFile, "Ostrich Run Record\n");
   fprintf(pFile, "%s", pBef);
   if(pObsGroup != NULL) pObsGroup->Write(pFile, WRITE_BNR, NULL);
   if(pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_BNR);
   algo->GetParamGroupPtr()->Write(pFile, WRITE_BNR);
   fprintf(pFile, "%s\n", pAft);
}/* end WriteBannerToFile() */

/******************************************************************************
WriteMultiObjRecord()

Write out multi-objective iteration result
******************************************************************************/
void WriteMultiObjRecord2(Algorithm *algo, int iter, ArchiveStruct * pArch, double dx)
{
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);
   pFile = fopen(fileName, "a");
   WriteMultiObjRecordToFile2(pFile, algo, iter, pArch, dx);
   fclose(pFile);

   sprintf(fileName, "OstNonDomSolutions%d.txt", 0);
   pFile = fopen(fileName, "w");
   WriteBannerToFile2(pFile, algo, "gen   ", "alg_conv_code");
   WriteMultiObjRecordToFile2(pFile, algo, iter, pArch, dx);
   fclose(pFile);

   WriteMultiObjRecordToFile2(stdout, algo, iter, pArch, dx);
}/* end WriteMultiObjRecord() */

/******************************************************************************
WriteMultiObjRecordToFile()

Write out multi-objective iteration result
******************************************************************************/
void WriteMultiObjRecordToFile2(FILE * pFile, Algorithm *algo, int iter, ArchiveStruct * pArch, double dx)
{
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   pGroup = algo->GetParamGroupPtr();

   fprintf(pFile, "\n");
   for(ArchiveStruct * pCur = pArch; pCur != NULL; pCur = pCur->pNext)
   {
      fprintf(pFile, "%-4d  ", iter);
      for(int i = 0; i < pCur->nF; i++)
      {
         fprintf(pFile, "%E  ", pCur->F[i]);
      }
      for(int i = 0; i < pCur->nX; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         fprintf(pFile, "%E  ", pParam->ConvertOutVal(pCur->X[i]));
      }
      fprintf(pFile, "%E\n", dx);
   }/* end for each non-dominated solution */
}/* end WriteMultiObjRecordToFile() */

/******************************************************************************
WriteRecord()

Write out iteration result
******************************************************************************/
void WriteRecord2(Algorithm *algo, int iter, double fx, double dx)
{
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);

   pFile = fopen(fileName, "a");
   WriteRecordToFile2(pFile, algo, iter, fx, dx);
   fclose(pFile);

   WriteRecordToFile2(stdout, algo, iter, fx, dx);
}/* end WriteRecord() */

/******************************************************************************
WriteRecordToFile()

Write out iteration result
******************************************************************************/
void WriteRecordToFile2(FILE * pFile, Algorithm *algo, int iter, double fx, double dx)
{
   ObservationGroup * pObsGroup;
   ResponseVarGroup * pRespVarGroup;
   ObjectiveFunction * pObjFunc;

   pObjFunc = algo->GetObjFuncPtr();
   pObsGroup = algo->GetObsGroupPtr();
   pRespVarGroup = NULL;
   if(pObjFunc != NULL)
      pRespVarGroup = (ResponseVarGroup *)(pObjFunc->GetResponseVarGroup());

   fprintf(pFile, "%-4d  %E  ", iter, fx);
   if(pObsGroup != NULL) pObsGroup->Write(pFile, WRITE_SCI, NULL);
   if(pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_SCI);
   algo->GetParamGroupPtr()->Write(pFile, WRITE_SCI);   
   fprintf(pFile, "%E\n", dx);
}/* end WriteRecordToFile() */

/******************************************************************************
WriteMultiObjOptimal()

Write out final set of dominated and non-dominated solutions for a multi-
objective application.
******************************************************************************/
void WriteMultiObjOptimal2(Algorithm *algo, ArchiveStruct * pNonDom, ArchiveStruct * pDom)
{
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);

   pFile = fopen(fileName, "a");
   WriteMultiObjOptimalToFile2(pFile, algo, pNonDom, pDom);
   //additional PATO output
   WriteCostToFile2(pFile, algo);
   WriteConstraintsToFile2(pFile, algo);
   WriteWellsToFile2(pFile, algo);   
   fclose(pFile);

   WriteMultiObjOptimalToFile2(stdout, algo, pNonDom, pDom);

   //additional PATO output
   WriteCostToFile2(stdout, algo);
   WriteConstraintsToFile2(stdout, algo);
   WriteWellsToFile2(stdout, algo);   
}/* end WriteMultiObjOptimal() */

/******************************************************************************
WriteMultiObjOptimalToFile()

Write out final set of dominated and non-dominated solutions for a multi-
objective application.
******************************************************************************/
void WriteMultiObjOptimalToFile2(FILE * pFile, Algorithm *algo, ArchiveStruct * pNonDom, ArchiveStruct * pDom)
{
   int nonDomCount = 0;
   int domCount = 0;

   ObservationGroup * pObsGroup;
   ResponseVarGroup * pRespVarGroup;
   ObjectiveFunction * pObjFunc;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   pObjFunc = algo->GetObjFuncPtr();
   pObsGroup = algo->GetObsGroupPtr();
   pGroup = algo->GetParamGroupPtr();
   pRespVarGroup = NULL;
   if(pObjFunc != NULL)
      pRespVarGroup = (ResponseVarGroup *)(pObjFunc->GetResponseVarGroup());

   fprintf(pFile, "\nNon-Dominated Solutions\n");

   //banner text
   if(pObsGroup != NULL) pObsGroup->Write(pFile, WRITE_BNR, NULL);
   if(pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_BNR);
   algo->GetParamGroupPtr()->Write(pFile, WRITE_BNR);   
   fprintf(pFile, "\n");

   //list of solutions
   for(ArchiveStruct * pCur = pNonDom; pCur != NULL; pCur = pCur->pNext)
   {
      for(int i = 0; i < pCur->nF; i++)
      {
         fprintf(pFile, "%E  ", pCur->F[i]);
      }
      for(int i = 0; i < pCur->nX; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         fprintf(pFile, "%E  ", pParam->ConvertOutVal(pCur->X[i]));
      }
      fprintf(pFile, "\n");
      nonDomCount++;
   }/* end for each non-dominated solution */
  
   fprintf(pFile, "\nDominated Solutions\n");
   //banner text
   if(pObsGroup != NULL) pObsGroup->Write(pFile, WRITE_BNR, NULL);
   if(pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_BNR);
   algo->GetParamGroupPtr()->Write(pFile, WRITE_BNR);   
   fprintf(pFile, "\n");

   //list of solutions
   for(ArchiveStruct * pCur = pDom; pCur != NULL; pCur = pCur->pNext)
   {
      for(int i = 0; i < pCur->nF; i++)
      {
         fprintf(pFile, "%E  ", pCur->F[i]);
      }
      for(int i = 0; i < pCur->nX; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         fprintf(pFile, "%E  ", pParam->ConvertOutVal(pCur->X[i]));
      }
      fprintf(pFile, "\n");
      domCount++;
   }/* end for each dominated solution */

   fprintf(pFile, "\nNumber of Non-Dominated Solutions : %d\n", nonDomCount);
   fprintf(pFile, "\nNumber of Dominated Solutions     : %d\n", domCount);
}/* end WriteMultiObjOptimalToFile() */

/******************************************************************************
WriteOptimal()

Write out optimal result
******************************************************************************/
void WriteOptimal2(Algorithm *algo, double fx)
{
   FILE * pFile;
   char fileName[DEF_STR_SZ];

   sprintf(fileName, "OstOutput%d.txt", 0);

   pFile = fopen(fileName, "a");
   WriteOptimalToFile2(pFile, algo, fx);
   //additional PATO output
   WriteCostToFile2(pFile, algo);
   WriteConstraintsToFile2(pFile, algo);
   WriteWellsToFile2(pFile, algo);   
   fclose(pFile);

   WriteOptimalToFile2(stdout, algo, fx);
   //additional PATO output
   WriteCostToFile2(stdout, algo);
   WriteConstraintsToFile2(stdout, algo);
   WriteWellsToFile2(stdout, algo);   
}/* end WriteOptimal() */

/******************************************************************************
WriteOptimalToFile()

Write out optimal result
******************************************************************************/
void WriteOptimalToFile2(FILE * pFile, Algorithm *algo, double fx) {
   fprintf(pFile, "\nOptimal Parameter Set\n");
   fprintf(pFile, "Objective Function : %E\n", fx);
   algo->GetParamGroupPtr()->Write(pFile, WRITE_OPT);
}/* end WriteOptimalToFile() */

/******************************************************************************
WriteOptimalToFileWithGroup()

Write out optimal result
******************************************************************************/
void WriteOptimalToFileWithGroup2(FILE* pFile, ParameterGroup* paramGroup, double fx) {

    fprintf(pFile, "\nOptimal Parameter Set\n");
    fprintf(pFile, "Objective Function : %E\n", fx);
    paramGroup->Write(pFile, WRITE_OPT);

}/* end WriteOptimalToFileWithGroup() */

/******************************************************************************
WriteAlgMetrics()

Write out algorithm metrics to stdout and to Ostrich output file.
******************************************************************************/
//void WriteAlgMetrics(AlgorithmABC * pAlg)
//{
//   FILE * pFile;
//   char fileName[DEF_STR_SZ];
//   int num_procs;
//   int id;
//
//   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
//   MPI_Comm_rank(MPI_COMM_WORLD, &id);
//
//   sprintf(fileName, "OstOutput%d.txt", id);
//
//   pFile = fopen(fileName, "a");
//   pAlg->WriteMetrics(pFile);
//   fclose(pFile);
//
//   pAlg->WriteMetrics(stdout);
//}/* end WriteMetrics() */

/******************************************************************************
WriteCostToFile()

Write out PATO cost breakdown.
******************************************************************************/
void WriteCostToFile2(FILE * pFile, Algorithm *algo)
{
   PATO * pPATO;
   if(algo->GetObjFuncId() != OBJ_FUNC_PATO){ return;}
   fprintf(pFile, "\nCost Breakdown\n");
   pPATO = (PATO *)(algo->GetObjFuncPtr());
   pPATO->WriteCost(pFile, WRITE_DEC);
}/* end WriteCostToFile() */

/******************************************************************************
WriteConstraintsToFile()

Write out PATO constraint information.
******************************************************************************/
void WriteConstraintsToFile2(FILE * pFile, Algorithm *algo)
{
   PATO * pPATO;
   GCOP * pGCOP;
   if(algo->GetObjFuncId() == OBJ_FUNC_PATO)
   { 
      fprintf(pFile, "\nSummary of Constraints\n");
      pPATO = (PATO *)(algo->GetObjFuncPtr());
      pPATO->WriteConstraints(pFile, WRITE_BNR);
      pPATO->WriteConstraints(pFile, WRITE_SCI);
   }
   else if (algo->GetObjFuncId() == OBJ_FUNC_GCOP)
   { 
      fprintf(pFile, "\nSummary of Constraints\n");
      pGCOP = (GCOP *)(algo->GetObjFuncPtr());
      pGCOP->WriteConstraints(pFile, WRITE_BNR);
      pGCOP->WriteConstraints(pFile, WRITE_SCI);
   }
   else { return;}
}/* end WriteConstraintsToFile() */

/******************************************************************************
WriteWellsToFile()

Write out PATO well information.
******************************************************************************/
void WriteWellsToFile2(FILE * pFile, Algorithm *algo)
{
   PATO * pPATO;
   if(algo->GetObjFuncId() != OBJ_FUNC_PATO){ return;}
   fprintf(pFile, "\nSummary of Optimal Wells\n");
   pPATO = (PATO *)(algo->GetObjFuncPtr());
   pPATO->WriteWells(pFile, WRITE_BNR);
   pPATO->WriteWells(pFile, WRITE_DEC);
}/* WriteWellsToFile() */

/******************************************************************************
WriteMelt()

Write out melting information. If input arg is less than or equal to zero, the 
banner is output, else the argument is treated as the melting count. The 'c' 
arg is used as an indicator of whether the melting operation increased or 
decreased the objective.
******************************************************************************/
void WriteMelt2(int count, int max, char c)
{
   FILE * pFile;
   char fileName[DEF_STR_SZ];
   int num_procs;
   int id;
   static bool CR = false; //true if last char was a newline

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   sprintf(fileName, "OstOutput%d.txt", id);

   pFile = fopen(fileName, "a");
   if(count == 0)
   {
      fprintf(pFile , "\nBeginning melting operation (requires %d evals)\n", max);
      fprintf(stdout, "\nBeginning melting operation (requires %d evals)\n", max);
      CR = true;
   }
   else if((count == -1) && (max == -1))
   {
      if(CR == false)
      {
         fprintf(pFile , "\n");
         fprintf(stdout, "\n");
      }
      fprintf(pFile , "Melting operation is complete\n\n");
      fprintf(stdout, "Melting operation is complete\n\n");
      CR = true;
   }
   else
   {
      fprintf(pFile , "%4d%c%c%c", count,c,c,c);
      fprintf(stdout, "%4d%c%c%c", count,c,c,c);

      CR = false;
      if((count % 10) == 0)
      {
         fprintf(pFile , "\n");
         fprintf(stdout, "\n");
         CR = true;
      }
   }
   fclose(pFile);
}/* end WriteMelt() */

/******************************************************************************
Write1dSearch()

Write out one-dimensional search information (to stdout only). If 'count' is:
   WRITE_GSECT : the Golden Section banner is output
   WRITE_BRENT : the Brent banner is output
   WRITE_SWTCH : indicates a switch from Brent to Golden Section
   WRITE_ENDED : end of 1d search
   >0          : treated as the search count.
******************************************************************************/
void Write1dSearch2(int count, int max)
{
   static bool NL = false;

   if(count == WRITE_GSECT)
   {
      fprintf(stdout, "\nBeginning Golden Section Search (requires %d evals)\n", max);
      NL = true;
   }
   else if(count == WRITE_BRENT)
   {
      fprintf(stdout, "\nBeginning Brent Search (max of %d evals)\n", max);
      NL = true;
   }
   else if(count == WRITE_SWTCH)
   {
      if(NL == false){fprintf(stdout, "\n");}
      fprintf(stdout, "\nGiving up on Brent method, switching to Golden Section\n");
      NL = true;
   }
   else if(count == WRITE_ENDED)
   {
      if(NL == false){fprintf(stdout, "\n");}
      fprintf(stdout, "Search operation is complete\n\n");
      NL = true;
   }
   else
   {
      fprintf(stdout, "%4d...", count);

      NL = false;
      if((count % 10) == 0)
      {    
         fprintf(stdout, "\n");
         NL = true;
      }
   }
}/* end Write1dSearch() */

/******************************************************************************
WriteInnerEval()

Write out inner loop information (to stdout only). If 'count' is:
   WRITE_GA    : the Genetic Algorithm banner is output
   WRITE_PSO   : the Particle Swarm Optimization banner is output
   WRITE_SMP   : the Sampling Algorithm banner is output
   WRITE_GRID  : the grid evaluation banner is output
   WRITE_SA    : the Simulated Annealing banner is output
   WRITE_LEV   : the Levenberg-Marquardt banner is output
   WRITE_DDS   : the Dynamically Dimensioned Search banner is output
   WRITE_USR   : the user defined evaluations banner is output
   WRITE_LHS   : the hypercube sampling banner is output
   WRITE_ENDED : end of inner evaluations
   WRITE_SCE   : the Shuffled Complex Evolution banner is output
   >0          : treated as the eval count.
******************************************************************************/
void WriteInnerEval2(int count, int max, char c)
{
   static bool LF = false;
   static int lfcount = 0;

   if(count == WRITE_BIS)
   {
      fprintf(stdout, "\nEvaluating inner bisections (requires at least %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_SMP)
   {
      fprintf(stdout, "\nEvaluating samples (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_JAC)
   {
      fprintf(stdout, "\nEvaluating global Jacobian (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_LHS)
   {
      fprintf(stdout, "\nEvaluating LHS samples (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_USR)
   {
      fprintf(stdout, "\nPerforming user-defined evaluations (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_GA)
   {
      fprintf(stdout, "\nEvaluating Population Fitness (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_PSO)
   {
      fprintf(stdout, "\nEvaluating Swarm (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_SCE)
   {
      fprintf(stdout, "\nEvaluating Complex (requires up to %d evals)\n", 3*max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_GLUE) //also used by Rejection Sampler
   {
      fprintf(stdout, "\nEvaluating samples (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_GRID)
   {
      fprintf(stdout, "\nEvaluating Mini Grid (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_SA)
   {
      fprintf(stdout, "\nPerforming Annealing Transitions (requires %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_LEV)
   {
      fprintf(stdout, "\nAdjusting Lambda Parameter (max of %d evals)\n", max);
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_ENDED)
   {
      if(LF == false){fprintf(stdout, "\n");}
      fprintf(stdout, "Operation is complete\n\n");
      LF = true;
      lfcount = 0;
   }
   else if(count == WRITE_DDS)
   {
      fprintf(stdout, "\nDDS is searching for a better solution.\n");
      LF = true;
      lfcount = 0;
   }
   else
   {
      fprintf(stdout, "%4d%c%c%c", count,c,c,c);
      LF = false;
      lfcount++;
      if(lfcount == 10)
      {    
         fprintf(stdout, "\n");
         LF = true;
         lfcount = 0;
      }
   }
}/* end WriteInnerEval() */

/******************************************************************************
WriteGrid()

Store parameter and objective function values to grid output file.
******************************************************************************/
void WriteGrid2(GridStruct * pGrid, int size)
{
   static int idx = 0;
   FILE * pFile;
   char name[DEF_STR_SZ];
   int id, i, j;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(name, "OstGrid%d.csv", id);

   if(idx == 0)
   {
      remove(name); 
      //write out banner.
      pFile = fopen(name, "w");
      fprintf(pFile,"Run,obj. function,");
      for(i = 0; i < pGrid->nprm; i++)
      {
         fprintf(pFile,"%s",GetParameterName(i));
         if(i < (pGrid->nprm - 1)) fprintf(pFile,",");
      }
      fprintf(pFile,"\n");
      fclose(pFile);
   }

   pFile = fopen(name, "a+");
   for(i = 0; i < size; i++)
   {
	   fprintf(pFile, "%d,%E,", idx++, pGrid->f[i]);
      for(j = 0; j < pGrid->nprm; j++)
      {
         fprintf(pFile,"%s",GetParameterValStr(j, pGrid->p[i][j]));
         if(j < (pGrid->nprm - 1)) fprintf(pFile,",");
      }
      fprintf(pFile,"\n");
   }
   fclose(pFile);
}/* end WriteGrid() */

/******************************************************************************
WritePreciseNumber()

Write a number (x) to the specified file (pFile) using the requested number
of digits of precision.
******************************************************************************/
void WritePreciseNumber2(FILE * pOut, double x)
{
   int precision = GetNumDigitsOfPrecision();
   switch(precision)
   {
      case (1) : fprintf(pOut, "%.1E", x); break;
      case (2) : fprintf(pOut, "%.2E", x); break;
      case (3) : fprintf(pOut, "%.3E", x); break;
      case (4) : fprintf(pOut, "%.4E", x); break;
      case (5) : fprintf(pOut, "%.5E", x); break;
      case (6) : fprintf(pOut, "%.6E", x); break;
      case (7) : fprintf(pOut, "%.7E", x); break;
      case (8) : fprintf(pOut, "%.8E", x); break;
      case (9) : fprintf(pOut, "%.9E", x); break;
      case (10) : fprintf(pOut, "%.10E", x); break;
      case (11) : fprintf(pOut, "%.11E", x); break;
      case (12) : fprintf(pOut, "%.12E", x); break;
      case (13) : fprintf(pOut, "%.13E", x); break;
      case (14) : fprintf(pOut, "%.14E", x); break;
      case (15) : fprintf(pOut, "%.15E", x); break;
      case (16) : fprintf(pOut, "%.16E", x); break;
      case (17) : fprintf(pOut, "%.17E", x); break;
      case (18) : fprintf(pOut, "%.18E", x); break;
      case (19) : fprintf(pOut, "%.19E", x); break;
      case (20) : fprintf(pOut, "%.20E", x); break;
      case (21) : fprintf(pOut, "%.21E", x); break;
      case (22) : fprintf(pOut, "%.22E", x); break;
      case (23) : fprintf(pOut, "%.23E", x); break;
      case (24) : fprintf(pOut, "%.24E", x); break;
      case (25) : fprintf(pOut, "%.25E", x); break;
      case (26) : fprintf(pOut, "%.26E", x); break;
      case (27) : fprintf(pOut, "%.27E", x); break;
      case (28) : fprintf(pOut, "%.28E", x); break;
      case (29) : fprintf(pOut, "%.29E", x); break;
      case (30) : fprintf(pOut, "%.30E", x); break;
      case (31) : fprintf(pOut, "%.31E", x); break;
      case (32) : fprintf(pOut, "%.32E", x); break;
      default  : fprintf(pOut, "%.6E", x); break;        
   }
}/* WritePreciseNumber() */
   
