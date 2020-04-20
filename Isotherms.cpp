/******************************************************************************
File     : Isotherms.cpp
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Isotherm classes are used to compute isotherms using the parameters and output 
concentrations specified in the IsothermIn.txt input file.

Version History
03-16-04    lsm   Created
03-24-04    lsm   Corrected naming of (1/n) Freundlich parameter
08-17-04    lsm   Added printing of RAM allocations
03-04-05    lsm   Added Polanyi-Partition isotherm
06-09-05    lsm   Added Langmuir-Partition isotherm
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
06-22-07    lsm      11. Dual Langmuir Isotherm
07-30-07    lsm      12. Orear Isotherm (for testing purposes only)
                     13. McCammon Isotherm (for testing purposes only)
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Isotherms.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObservationGroup.h"
#include "Observation.h"

#include "IsoParse.h"

/******************************************************************************
LinearIsotherm::Initialize()

Initialize parameters and output arrays, using input string.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LinearIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Kd = 0.00;
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginLinearIsotherm")         == NULL){ strcat(msg, "BeginLinearIsotherm\n");}
   if(strstr(pStr, "EndLinearIsotherm")           == NULL){ strcat(msg, "EndLinearIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Linear Isotherm section
   pTmp = strstr(pStr,"BeginLinearIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndLinearIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "Kd") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Kd);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end LinearIsotherm::Initialize() */

/******************************************************************************
LinearIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LinearIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Kd = pgroup->GetParamPtr("KdVal")->GetTransformedVal();
   return true;
}/* end LinearIsotherm::Initialize() */

/******************************************************************************
LinearIsotherm::Destroy()

Free up memory.
******************************************************************************/
void LinearIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end LinearIsotherm::Destroy() */

/******************************************************************************
LinearIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void LinearIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Linear\n");
   fprintf(pFile, "Kd %E\n", m_Kd);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end LinearIsotherm::Compute() */

/******************************************************************************
LangmuirIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LangmuirIsotherm::Initialize(char * pStr)
{
   bool bLumpedQ0 = false;
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   double Q0_b = 0.00;
   m_Q0 = 0.00;
   m_b = 0.00;
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginLangmuirIsotherm")         == NULL){ strcat(msg, "BeginLangmuirIsotherm\n");}
   if(strstr(pStr, "EndLangmuirIsotherm")           == NULL){ strcat(msg, "EndLangmuirIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Langmuir Isotherm section
   pTmp = strstr(pStr,"BeginLangmuirIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndLangmuirIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "b*Q0") != NULL)
      {
         sscanf(pLine, "%s %lf", pVar, &Q0_b);
         bLumpedQ0 = true;
      }
      else if(strstr(pLine, "b") != NULL) {sscanf(pLine, "%s %lf", pVar, &m_b);}
      else if(strstr(pLine, "Q0") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Q0);}
   }

   if(bLumpedQ0 == true) { m_Q0 = (Q0_b/m_b);}

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end LangmuirIsotherm::Initialize() */

/******************************************************************************
LangmuirIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LangmuirIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   return true;
}/* end LangmuirIsotherm::Initialize() */

/******************************************************************************
LangmuirIsotherm::Destroy()

Free up memory.
******************************************************************************/
void LangmuirIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end LangmuirIsotherm::Destroy() */

/******************************************************************************
LangmuirIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void LangmuirIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Langmuir\n");
   fprintf(pFile, "b*Q0 %E\n", m_Q0*m_b);
   fprintf(pFile, "b    %E\n", m_b);
   fprintf(pFile, "Q0   %E\n", m_Q0);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end LangmuirIsotherm::Compute() */

/******************************************************************************
DualLangmuirIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool DualLangmuirIsotherm::Initialize(char * pStr)
{
   bool bLumpedQ0 = false;
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   double Q01_b1, Q02_b2;
   Q01_b1 = Q02_b2 = 0.00;
   m_Q01 = 0.00;
   m_b1  = 0.00;
   m_Q02 = 0.00;
   m_b2  = 0.00;

   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginDualLangmuirIsotherm")         == NULL){ strcat(msg, "BeginDualLangmuirIsotherm\n");}
   if(strstr(pStr, "EndDualLangmuirIsotherm")           == NULL){ strcat(msg, "EndDualLangmuirIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Langmuir Isotherm section
   pTmp = strstr(pStr,"BeginDualLangmuirIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndDualLangmuirIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "b1*Q01") != NULL)
      {
         sscanf(pLine, "%s %lf", pVar, &Q01_b1);
         bLumpedQ0 = true;
      }
      else if(strstr(pLine, "b1") != NULL) {sscanf(pLine, "%s %lf", pVar, &m_b1);}
      else if(strstr(pLine, "Q01") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Q01);}
      else if(strstr(pLine, "b2*Q02") != NULL)
      {
         sscanf(pLine, "%s %lf", pVar, &Q02_b2);
         bLumpedQ0 = true;
      }
      else if(strstr(pLine, "b2") != NULL) {sscanf(pLine, "%s %lf", pVar, &m_b2);}
      else if(strstr(pLine, "Q02") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Q02);}
   }

   if(bLumpedQ0 == true) 
   {
      m_Q01 = (Q01_b1/m_b1);
      m_Q02 = (Q02_b2/m_b2);
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end DualLangmuirIsotherm::Initialize() */

/******************************************************************************
DualLangmuirIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool DualLangmuirIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q01 = pgroup->GetParamPtr("Q01Val")->GetTransformedVal();
   m_b1 = pgroup->GetParamPtr("b1Val")->GetTransformedVal();
   m_Q02 = pgroup->GetParamPtr("Q02Val")->GetTransformedVal();
   m_b2 = pgroup->GetParamPtr("b2Val")->GetTransformedVal();
   return true;
}/* end DualLangmuirIsotherm::Initialize() */

/******************************************************************************
DualLangmuirIsotherm::Destroy()

Free up memory.
******************************************************************************/
void DualLangmuirIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end DualLangmuirIsotherm::Destroy() */

/******************************************************************************
DualLangmuirIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void DualLangmuirIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Dual-Langmuir\n");
   fprintf(pFile, "b1*Q01 %E\n", m_Q01*m_b1);
   fprintf(pFile, "b1    %E\n", m_b1);
   fprintf(pFile, "Q01   %E\n", m_Q01);
   fprintf(pFile, "b2*Q02 %E\n", m_Q02*m_b2);
   fprintf(pFile, "b2    %E\n", m_b2);
   fprintf(pFile, "Q02   %E\n", m_Q02);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end DualLangmuirIsotherm::Compute() */

/******************************************************************************
FreundlichIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool FreundlichIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Kf = 0.00;
   m_Nf = 0.00;
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginFreundlichIsotherm")     == NULL){ strcat(msg, "BeginFreundlichIsotherm\n");}
   if(strstr(pStr, "EndFreundlichIsotherm")       == NULL){ strcat(msg, "EndFreundlichIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Linear Isotherm section
   pTmp = strstr(pStr,"BeginFreundlichIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndFreundlichIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "Kf") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Kf);}
      else if(strstr(pLine, "(1/n)") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Nf);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end FreundlichIsotherm::Initialize() */

/******************************************************************************
FreundlichIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool FreundlichIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Kf = pgroup->GetParamPtr("KfVal")->GetTransformedVal();
   m_Nf = pgroup->GetParamPtr("(1/n)Val")->GetTransformedVal();
   return true;
}/* end FreundlichIsotherm::Initialize() */

/******************************************************************************
FreundlichIsotherm::Destroy()

Free up memory.
******************************************************************************/
void FreundlichIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end FreundlichIsotherm::Destroy() */

/******************************************************************************
FreundlichIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void FreundlichIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Freundlich\n");
   fprintf(pFile, "Kf  %E\n", m_Kf);
   fprintf(pFile, "(1/n)  %E\n", m_Nf);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end FreundlichIsotherm::Compute() */

/******************************************************************************
FreundlichPartitionIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool FreundlichPartitionIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Kf = 0.00;
   m_Nf = 0.00;
   m_Kp = 0.00;
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginFreundlich-PartitionIsotherm")     == NULL){ strcat(msg, "BeginFreundlich-PartitionIsotherm\n");}
   if(strstr(pStr, "EndFreundlich-PartitionIsotherm")       == NULL){ strcat(msg, "EndFreundlich-PartitionIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Linear Isotherm section
   pTmp = strstr(pStr,"BeginFreundlich-PartitionIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndFreundlich-PartitionIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "Kf") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Kf);}
      else if(strstr(pLine, "Kp") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Kp);}
      else if(strstr(pLine, "(1/n)") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Nf);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end FreundlichPartitionIsotherm::Initialize() */

/******************************************************************************
FreundlichPartitionIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool FreundlichPartitionIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Kf = pgroup->GetParamPtr("KfVal")->GetTransformedVal();
   m_Kp = pgroup->GetParamPtr("KpVal")->GetTransformedVal();
   m_Nf = pgroup->GetParamPtr("(1/n)Val")->GetTransformedVal();
   return true;
}/* end FreundlichPartitionIsotherm::Initialize() */

/******************************************************************************
FreundlichPartitionIsotherm::Destroy()

Free up memory.
******************************************************************************/
void FreundlichPartitionIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end FreundlichPartitionIsotherm::Destroy() */

/******************************************************************************
FreundlichPartitionIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void FreundlichPartitionIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Freundlich-Partition\n");
   fprintf(pFile, "Kp  %E\n", m_Kp);
   fprintf(pFile, "Kf  %E\n", m_Kf);
   fprintf(pFile, "(1/n)  %E\n", m_Nf);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end FreundlichPartitionIsotherm::Compute() */

/******************************************************************************
PolanyiIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool PolanyiIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Q0 = 0.00;
   m_b = 0.00;
   m_a = 0.00;
   m_Sw = 0.00;

   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginPolanyiIsotherm") == NULL){ strcat(msg, "BeginPolanyiIsotherm\n");}
   if(strstr(pStr, "EndPolanyiIsotherm")   == NULL){ strcat(msg, "EndPolanyiIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Isotherm section
   pTmp = strstr(pStr,"BeginPolanyiIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndPolanyiIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "Q0") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Q0);}
      else if(strstr(pLine, "a") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_a);}
      else if(strstr(pLine, "b") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_b);}
      else if(strstr(pLine, "Sw") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Sw);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end PolanyiIsotherm::Initialize() */

/******************************************************************************
PolanyiIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool PolanyiIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_a = pgroup->GetParamPtr("aVal")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   return true;
}/* end PolanyiIsotherm::Initialize() */

/******************************************************************************
PolanyiIsotherm::Destroy()

Free up memory.
******************************************************************************/
void PolanyiIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end PolanyiIsotherm::Destroy() */

/******************************************************************************
PolanyiIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void PolanyiIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Polanyi\n");
   fprintf(pFile, "Sw  %E\n", m_Sw);
   fprintf(pFile, "Q0  %E\n", m_Q0);
   fprintf(pFile, "a   %E\n", m_a);
   fprintf(pFile, "b   %E\n", m_b);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end PolanyiIsotherm::Compute() */

/******************************************************************************
PolanyiPartitionIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool PolanyiPartitionIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Kp = 0.00;
   m_Q0 = 0.00;
   m_b = 0.00;
   m_a = 0.00;
   m_Sw = 0.00;

   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginPolanyi-PartitionIsotherm") == NULL){ strcat(msg, "BeginPolanyi-PartitionIsotherm\n");}
   if(strstr(pStr, "EndPolanyi-PartitionIsotherm")   == NULL){ strcat(msg, "EndPolanyi-PartitionIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Isotherm section
   pTmp = strstr(pStr,"BeginPolanyi-PartitionIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndPolanyi-PartitionIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "Kp") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Kp);}
      else if(strstr(pLine, "Q0") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Q0);}
      else if(strstr(pLine, "a") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_a);}
      else if(strstr(pLine, "b") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_b);}
      else if(strstr(pLine, "Sw") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Sw);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end PolanyiPartitionIsotherm::Initialize() */

/******************************************************************************
PolanyiPartitionIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool PolanyiPartitionIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_Kp = pgroup->GetParamPtr("KpVal")->GetTransformedVal();
   m_a = pgroup->GetParamPtr("aVal")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   return true;
}/* end PolanyiPartitionIsotherm::Initialize() */

/******************************************************************************
PolanyiPartitionIsotherm::Destroy()

Free up memory.
******************************************************************************/
void PolanyiPartitionIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end PolanyiPartitionIsotherm::Destroy() */

/******************************************************************************
PolanyiPartitionIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void PolanyiPartitionIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Polanyi-Partition\n");
   fprintf(pFile, "Sw  %E\n", m_Sw);
   fprintf(pFile, "Kp  %E\n", m_Kp);
   fprintf(pFile, "Q0  %E\n", m_Q0);
   fprintf(pFile, "a   %E\n", m_a);
   fprintf(pFile, "b   %E\n", m_b);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end PolanyiPartitionIsotherm::Compute() */

/******************************************************************************
LangmuirPartitionIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LangmuirPartitionIsotherm::Initialize(char * pStr)
{
   int i;
   bool bLumpedQ0 = false;
   char * pTmp;
   
   const char * begin = "BeginLangmuir-PartitionIsotherm";
   const char * end   = "EndLangmuir-PartitionIsotherm";

   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Kp = 0.00;
   double Q0_b = 0.00;
   m_Q0 = 0.00;
   m_b = 0.00;
   
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, begin) == NULL){ strcat(msg, begin); strcat(msg, "\n");}
   if(strstr(pStr, end)   == NULL){ strcat(msg, end); strcat(msg, "\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Isotherm section
   pTmp = strstr(pStr,begin);
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, end) != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "Kp") != NULL){sscanf(pLine, "%s %lf", pVar, &m_Kp);}
      else if(strstr(pLine, "b*Q0") != NULL)
      { 
         sscanf(pLine, "%s %lf", pVar, &Q0_b);
         bLumpedQ0 = true;
      }
      else if(strstr(pLine, "b") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_b);}
      else if(strstr(pLine, "Q0") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Q0);}
   }

   if(bLumpedQ0 == true) { m_Q0 = (Q0_b/m_b);}

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end LangmuirPartitionIsotherm::Initialize() */

/******************************************************************************
LangmuirPartitionIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LangmuirPartitionIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_Kp = pgroup->GetParamPtr("KpVal")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   return true;
}/* end LangmuirPartitionIsotherm::Initialize() */

/******************************************************************************
LangmuirPartitionIsotherm::Destroy()

Free up memory.
******************************************************************************/
void LangmuirPartitionIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end LangmuirPartitionIsotherm::Destroy() */

/******************************************************************************
LangmuirPartitionIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void LangmuirPartitionIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Langmuir-Partition\n");
   fprintf(pFile, "Kp    %E\n", m_Kp);
   fprintf(pFile, "b*Q0  %E\n", m_Q0*m_b);
   fprintf(pFile, "b     %E\n", m_b);
   fprintf(pFile, "Q0    %E\n", m_Q0);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end LangmuirPartitionIsotherm::Compute() */

/******************************************************************************
BET_Isotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool BET_Isotherm::Initialize(char * pStr)
{
   int i;
   bool bLumpedQ0 = false;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   double Q0_b = 0.00;
   m_Q0 = 0.00;
   m_b = 0.00;
   m_Sw = 0.00;

   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginBET_Isotherm")     == NULL){ strcat(msg, "BeginBET_Isotherm\n");}
   if(strstr(pStr, "EndBET_Isotherm")       == NULL){ strcat(msg, "EndBET_Isotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Isotherm section
   pTmp = strstr(pStr,"BeginBET_Isotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndBET_Isotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "b*Q0") != NULL)
      {
         sscanf(pLine, "%s %lf", pVar, &Q0_b);
         bLumpedQ0 = true;
      }
      else if(strstr(pLine, "b") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_b);}
      else if(strstr(pLine, "Sw") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Sw);}
      else if(strstr(pLine, "Q0") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Q0);}
   }

   if(bLumpedQ0 == true) { m_Q0 = (Q0_b/m_b);}

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end BET_Isotherm::Initialize() */

/******************************************************************************
BET_Isotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool BET_Isotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   return true;
}/* end BET_Isotherm::Initialize() */

/******************************************************************************
BET_Isotherm::Destroy()

Free up memory.
******************************************************************************/
void BET_Isotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end BET_Isotherm::Destroy() */

/******************************************************************************
BET_Isotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void BET_Isotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType BET\n");
   fprintf(pFile, "Sw   %E\n", m_Sw);
   fprintf(pFile, "b*Q0 %E\n", m_Q0*m_b);
   fprintf(pFile, "b    %E\n", m_b);
   fprintf(pFile, "Q0   %E\n", m_Q0);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end BET_Isotherm::Compute() */

/******************************************************************************
TothIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool TothIsotherm::Initialize(char * pStr)
{
   int i;
   bool bLumpedQ0 = false;
   char * pTmp;
   
   const char * begin = "BeginTothIsotherm";
   const char * end   = "EndTothIsotherm";

   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   double Q0_b = 0.00;
   m_Q0 = 0.00;
   m_nT = 0.00;
   m_b = 0.00;
   
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, begin) == NULL){ strcat(msg, begin); strcat(msg, "\n");}
   if(strstr(pStr, end)   == NULL){ strcat(msg, end); strcat(msg, "\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Isotherm section
   pTmp = strstr(pStr,begin);
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, end) != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "n") != NULL){sscanf(pLine, "%s %lf", pVar, &m_nT);}
      else if(strstr(pLine, "b*Q0") != NULL)
      { 
         sscanf(pLine, "%s %lf", pVar, &Q0_b);
         bLumpedQ0 = true;
      }
      else if(strstr(pLine, "b") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_b);}
      else if(strstr(pLine, "Q0") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Q0);}
   }

   if(bLumpedQ0 == true){m_Q0 = (Q0_b/m_b);}

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end TothIsotherm::Initialize() */

/******************************************************************************
TothIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool TothIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   m_nT = pgroup->GetParamPtr("nVal")->GetTransformedVal();
   return true;
}/* end TothIsotherm::Initialize() */

/******************************************************************************
TothIsotherm::Destroy()

Free up memory.
******************************************************************************/
void TothIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end TothIsotherm::Destroy() */

/******************************************************************************
TothIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void TothIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Toth\n");
   fprintf(pFile, "b*Q0  %E\n", m_Q0*m_b);
   fprintf(pFile, "b     %E\n", m_b);
   fprintf(pFile, "n     %E\n", m_nT);
   fprintf(pFile, "Q0    %E\n", m_Q0);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end TothIsotherm::Compute() */

/******************************************************************************
LangmuirFreundlichIsotherm::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LangmuirFreundlichIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   
   const char * begin = "BeginLangmuir-FreundlichIsotherm";
   const char * end   = "EndLangmuir-FreundlichIsotherm";

   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_Q0 = 0.00;
   m_nG = 0.00;
   m_b = 0.00;
   
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, begin) == NULL){ strcat(msg, begin); strcat(msg, "\n");}
   if(strstr(pStr, end)   == NULL){ strcat(msg, end); strcat(msg, "\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Isotherm section
   pTmp = strstr(pStr,begin);
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, end) != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "(1/n)") != NULL){sscanf(pLine, "%s %lf", pVar, &m_nG);}
      else if(strstr(pLine, "Q0") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_Q0);}
      else if(strstr(pLine, "b") != NULL){ sscanf(pLine, "%s %lf", pVar, &m_b);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end LangmuirFreundlichIsotherm::Initialize() */

/******************************************************************************
LangmuirFreundlichIsotherm::Initialize()

Initialize parameter values using ParameterGroup.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool LangmuirFreundlichIsotherm::Initialize(ParameterGroup * pgroup)
{
   if(pgroup == NULL) return false;
   m_Q0 = pgroup->GetParamPtr("Q0Val")->GetTransformedVal();
   m_b = pgroup->GetParamPtr("bVal")->GetTransformedVal();
   m_nG = pgroup->GetParamPtr("(1/n)Val")->GetTransformedVal();
   return true;
}/* end LangmuirFreundlichIsotherm::Initialize() */

/******************************************************************************
LangmuirFreundlichIsotherm::Destroy()

Free up memory.
******************************************************************************/
void LangmuirFreundlichIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end LangmuirFreundlichIsotherm::Destroy() */

/******************************************************************************
LangmuirFreundlichIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void LangmuirFreundlichIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Langmuir-Freundlich\n");
   fprintf(pFile, "Q0    %E\n", m_Q0);
   fprintf(pFile, "b     %E\n", m_b);
   fprintf(pFile, "(1/n) %E\n", m_nG);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end LangmuirFreundlichIsotherm::Compute() */

/******************************************************************************
OrearIsotherm::Initialize()

Initialize parameters and output arrays, using input string.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool OrearIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_a = m_b = 0.00;
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginOrearIsotherm")         == NULL){ strcat(msg, "BeginOrearIsotherm\n");}
   if(strstr(pStr, "EndOrearIsotherm")           == NULL){ strcat(msg, "EndOrearIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Orear Isotherm section
   pTmp = strstr(pStr,"BeginOrearIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndOrearIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "a") != NULL){sscanf(pLine, "%s %lf", pVar, &m_a);}
      if(strstr(pLine, "b") != NULL){sscanf(pLine, "%s %lf", pVar, &m_b);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end OrearIsotherm::Initialize() */

/******************************************************************************
OrearIsotherm::Destroy()

Free up memory.
******************************************************************************/
void OrearIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end OrearIsotherm::Destroy() */

/******************************************************************************
OrearIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void OrearIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType Orear\n");
   fprintf(pFile, "a %E\n", m_a);
   fprintf(pFile, "b %E\n", m_b);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end OrearIsotherm::Compute() */

/******************************************************************************
McCammonIsotherm::Initialize()

Initialize parameters and output arrays, using input string.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool McCammonIsotherm::Initialize(char * pStr)
{
   int i;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char msg[DEF_STR_SZ];
   char * pLine;
   strcpy(m_OutFile, ISO_OUT_FILE);
   m_A = m_B = m_C = m_E = 0.00;
   m_NumOut = 0;
   m_pC = NULL;
   m_pq = NULL;

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginMcCammonIsotherm")         == NULL){ strcat(msg, "BeginMcCammonIsotherm\n");}
   if(strstr(pStr, "EndMcCammonIsotherm")           == NULL){ strcat(msg, "EndMcCammonIsotherm\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Orear Isotherm section
   pTmp = strstr(pStr,"BeginMcCammonIsotherm");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndMcCammonIsotherm") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "A") != NULL){sscanf(pLine, "%s %lf", pVar, &m_A);}
      if(strstr(pLine, "B") != NULL){sscanf(pLine, "%s %lf", pVar, &m_B);}
      if(strstr(pLine, "C") != NULL){sscanf(pLine, "%s %lf", pVar, &m_C);}
      if(strstr(pLine, "_E_") != NULL){sscanf(pLine, "%s %lf", pVar, &m_E);}
   }

   //parse the Concentrations section
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      m_NumOut++;
   }
   m_NumOut--;

   NEW_PRINT("double", m_NumOut);
   m_pC = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf", pVar, &(m_pC[i]));
   }
   return true;
}/* end McCammonIsotherm::Initialize() */

/******************************************************************************
McCammonIsotherm::Destroy()

Free up memory.
******************************************************************************/
void McCammonIsotherm::Destroy(void)
{
   delete [] m_pC;
   delete [] m_pq;
}/* end McCammonIsotherm::Destroy() */

/******************************************************************************
McCammonIsotherm::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void McCammonIsotherm::Compute(void)
{
   FILE * pFile;
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   pFile = fopen(m_OutFile, "w");
   fprintf(pFile, "Isotherm Output File\n\n");
   fprintf(pFile, "IsothermType McCammon\n");
   fprintf(pFile, "A %E\n", m_A);
   fprintf(pFile, "B %E\n", m_B);
   fprintf(pFile, "C %E\n", m_C);
   fprintf(pFile, "D %E (fixed)\n", m_D);
   fprintf(pFile, "E %E\n", m_E);
   fprintf(pFile, "F %E (fixed)\n", m_F);
   fprintf(pFile, "G %E (fixed)\n", -1.00);
   fprintf(pFile, "NumObs %d\n\n", m_NumOut);
   fprintf(pFile, "i     Concentration  q\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%-4d  %.16E  %.16E\n",i,m_pC[i],m_pq[i]);
   }   
   fclose(pFile);
}/* end McCammonIsotherm::Compute() */

/******************************************************************************
McCammonIsotherm::q()

Adapted from APPENDIX of McCammon, 1973, Mathematical Geology, vol. 5, no. 4, 
pg. 375.
******************************************************************************/
double McCammonIsotherm::q(double c)
{
   double YS;  //q
   double YDS; //dqdc
   double XS=c;
   double P1=m_A;
   double P2=m_B;
   double P3=m_C;
   double P4=m_D;
   double P5=m_E;
   double P6=m_F;
   double P7=-1.00;
   double P8=NEARLY_HUGE;
   double CE=2*P3*XS+P5;
   double C=P3*XS*XS+P5*XS+P6;
   double B=P2*XS+P4;
   double D=B*B-4*P1*C;
   double DSQ;

   //Test for zero discriminant
   if(D >= 0) goto ten;
   YDS=YS=P8;
   goto sixty;
ten:
   //Test for non-zero A
   if(P1 != 0) goto thirty;
   //Test for zero A and non-zero (Bx+D)
   if((P1 == 0.00) && (B != 0.00)) goto twenty;
   YDS=YS=P8;
   goto sixty;
twenty:
   YS=-C/B;
   YDS=(-CE*B+P2*C)/(B*B);
   goto sixty;
thirty:
   DSQ=sqrt(D);
   YS=(-B+P7*DSQ)/(2*P1);
   YDS=(-P2+P7*(P2*B-2*P1*CE)/DSQ)/(2*P1);
sixty:
   return YS;
}/* end McCammonIsotherm::q() */

/******************************************************************************
McCammonIsotherm::dqdc()

Adapted from APPENDIX of McCammon, 1973, Mathematical Geology, vol. 5, no. 4, 
pg. 375.
******************************************************************************/
double McCammonIsotherm::dqdc(double c)
{
   double YS;  //q
   double YDS; //dqdc
   double XS=c;
   double P1=m_A;
   double P2=m_B;
   double P3=m_C;
   double P4=m_D;
   double P5=m_E;
   double P6=m_F;
   double P7=-1.00;
   double P8=NEARLY_HUGE;
   double CE=2*P3*XS+P5;
   double C=P3*XS*XS+P5*XS+P6;
   double B=P2*XS+P4;
   double D=B*B-4*P1*C;
   double DSQ;

   //Test for zero discriminant
   if(D >= 0) goto ten;
   YDS=YS=P8;
   goto sixty;
ten:
   //Test for non-zero A
   if(P1 != 0) goto thirty;
   //Test for zero A and non-zero (Bx+D)
   if((P1 == 0.00) && (B != 0.00)) goto twenty;
   YDS=YS=P8;
   goto sixty;
twenty:
   YS=-C/B;
   YDS=(-CE*B+P2*C)/(B*B);
   goto sixty;
thirty:
   DSQ=sqrt(D);
   YS=(-B+P7*DSQ)/(2*P1);
   YDS=(-P2+P7*(P2*B-2*P1*CE)/DSQ)/(2*P1);
sixty:
   return YDS;
}/* end McCammonIsotherm::dqdc() */

/******************************************************************************
LinearIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void LinearIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end LinearIsotherm::Compute() */

/******************************************************************************
LangmuirIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void LangmuirIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end LangmuirIsotherm::Compute() */

/******************************************************************************
LangmuirPartitionIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void LangmuirPartitionIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end LangmuirPartitionIsotherm::Compute() */

/******************************************************************************
DualLangmuirIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void DualLangmuirIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end DualLangmuirIsotherm::Compute() */

/******************************************************************************
FreundlichIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void FreundlichIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end FreundlichIsotherm::Compute() */

/******************************************************************************
FreundlichPartitionIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void FreundlichPartitionIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end FreundlichPartitionIsotherm::Compute() */

/******************************************************************************
LangmuirFreundlichIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void LangmuirFreundlichIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end LangmuirFreundlichIsotherm::Compute() */

/******************************************************************************
PolanyiIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void PolanyiIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end PolanyiIsotherm::Compute() */

/******************************************************************************
PolanyiPartitionIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void PolanyiPartitionIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end PolanyiPartitionIsotherm::Compute() */

/******************************************************************************
TothIsotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void TothIsotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end TothIsotherm::Compute() */

/******************************************************************************
BET_Isotherm::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void BET_Isotherm::Compute(ObservationGroup * ogroup)
{
   int i;
   for(i = 0; i < m_NumOut; i++){ m_pq[i] = q(m_pC[i]);}
   for(i = 0; i < m_NumOut; i++){ ogroup->GetObsPtr(i)->SetComputedVal(m_pq[i]); }
}/* end BET_Isotherm::Compute() */



