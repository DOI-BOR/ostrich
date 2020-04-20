/******************************************************************************
File      : TiedRespVar.cpp
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

Contains defintions for 'tied' response variables. Tied response variables 
are computed as functions of one or more response variables, which are read 
from model input and/or output files. The ABC for response variable 
encapsulates the interface used by other Ostrich modules, allowing various 
specific tied response variable relationships (linear, exponential, etc.) to 
be implemented as needed with minimal code change (just need to add the 
specific tied response variable class and some additional input file parsing).

These specific tied-response variable classes are supported:
TiedRespVarLin1  : linear function of one response variable
TiedRespVarLin2  : linear function of two response variables
TiedRespVarWsum  : weighted sum of one or more response variables

Version History
01-10-05    lsm   Created
01-20-05    lsm   added support for weighted sum tied resp. vars.
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "TiedRespVar.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
TiedRespVarLin1::Destroy()
******************************************************************************/
void TiedRespVarLin1::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedRespVarLin1)
******************************************************************************/
TiedRespVarLin1::TiedRespVarLin1(void)
{
   m_pName = NULL;
   m_pTie  = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedRespVarLin1)
******************************************************************************/
TiedRespVarLin1::TiedRespVarLin1(IroncladString name, RespVarABC * p1, 
                             UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie = p1;  

   //parse the config string to determine value for C0 and C1
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedRespVarLin1()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetCurrentVal()

Computes the current value of the tied resp. var.
******************************************************************************/
double TiedRespVarLin1::GetCurrentVal(void)
{
   double y, x;

   x = m_pTie->GetCurrentVal();
   y = m_C1*x + m_C0;
   return y;   
} /* end GetCurrentVal() */

/******************************************************************************
GetInitialVal()

Computes the initial value of the tied resp. var.
******************************************************************************/
double TiedRespVarLin1::GetInitialVal(void)
{
   double y, x;

   x = m_pTie->GetInitialVal();
   y = m_C1*x + m_C0;
   return y;   
} /* end GetInitialVal() */

/******************************************************************************
TiedRespVarLin1::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedRespVarLin1::Write(FILE * pFile, int type)
{
   double val;

   val = GetCurrentVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Resp. Var. = %s\n", m_pTie->GetName());
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "Value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end TiedRespVarLin1::WriteToFile() */

/******************************************************************************
TiedRespVarLin2::Destroy()
******************************************************************************/
void TiedRespVarLin2::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedRespVarLin2)
******************************************************************************/
TiedRespVarLin2::TiedRespVarLin2(void)
{
   m_pName = NULL;
   m_pTie1  = NULL;
   m_pTie2  = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedRespVarLin2)
******************************************************************************/
TiedRespVarLin2::TiedRespVarLin2(IroncladString name, RespVarABC * p1, 
                             RespVarABC * p2, UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie1 = p1;
   m_pTie2 = p2;

   //parse the config string to determin valuese for C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedRespVarLin2()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedRespVarLin2()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedRespVarLin2()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetCurrentVal()

Computes the current value of the tied resp. var.
******************************************************************************/
double TiedRespVarLin2::GetCurrentVal(void)
{
   double y, x1, x2;

   x1 = m_pTie1->GetCurrentVal();
   x2 = m_pTie2->GetCurrentVal();
   y = m_C3*x1*x2 + m_C2*x2 + m_C1*x1 + m_C0;

   return y;   
} /* end GetCurrentVal() */

/******************************************************************************
GetInitialVal()

Computes the initial value of the tied resp. var.
******************************************************************************/
double TiedRespVarLin2::GetInitialVal(void)
{
   double y, x1, x2;

   x1 = m_pTie1->GetInitialVal();
   x2 = m_pTie2->GetInitialVal();
   y = m_C3*x1*x2 + m_C2*x2 + m_C1*x1 + m_C0;

   return y;   
} /* end GetInitialVal() */

/******************************************************************************
TiedRespVarLin2::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedRespVarLin2::Write(FILE * pFile, int type)
{
   double val;

   val = GetCurrentVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Resp. Var. #1 = %s\n", m_pTie1->GetName());
      fprintf(pFile, "Tied Resp. Var. #2 = %s\n", m_pTie2->GetName());
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
      fprintf(pFile, "C3 = %lf\n", m_C3);
      fprintf(pFile, "value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end TiedRespVarLin2::WriteToFile() */

/******************************************************************************
TiedRespVarWsum::Destroy()
******************************************************************************/
void TiedRespVarWsum::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pList;
   delete [] m_pWgt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedRespVarWsum)
******************************************************************************/
TiedRespVarWsum::TiedRespVarWsum(void)
{
   m_pName = NULL;
   m_pWgt  = NULL;
   m_pList = NULL;
   m_Num = 0;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedRespVarWsum)
******************************************************************************/
TiedRespVarWsum::TiedRespVarWsum(IroncladString name, RespVarABC ** pList, int nrv,
                             UnmoveableString configStr)
{
   int i, j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   NEW_PRINT("RespVarABC *", nrv);
   m_pList = new RespVarABC *[nrv];
   MEM_CHECK(m_pList);
   for(i = 0; i < nrv; i++) m_pList[i] = pList[i];

   m_Num = nrv;

   NEW_PRINT("double", nrv);
   m_pWgt = new double[nrv];
   MEM_CHECK(m_pWgt);
   for(i = 0; i < nrv; i++) m_pWgt[i] = 0.00;

   //parse the config string to determine value for weights
   pTok = configStr;
   for(i = 0; i < nrv; i++)
   {
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, i, nrv, "TiedRespVarWsum()");
      m_pWgt[i] = atof(tmpStr);
      pTok += j;
   }
      
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetCurrentVal()

Computes the current value of the tied resp. var.
******************************************************************************/
double TiedRespVarWsum::GetCurrentVal(void)
{
   int i;
   double wsum, w, x;

   wsum = 0.00;

   for(i = 0; i < m_Num; i++)
   {
      x = m_pList[i]->GetCurrentVal();
      w = m_pWgt[i];
      wsum += w*x;
   }
   return wsum;   
} /* end GetCurrentVal() */

/******************************************************************************
GetInitialVal()

Computes the initial value of the tied resp. var.
******************************************************************************/
double TiedRespVarWsum::GetInitialVal(void)
{
   int i;
   double wsum, w, x;

   wsum = 0.00;

   for(i = 0; i < m_Num; i++)
   {
      x = m_pList[i]->GetInitialVal();
      w = m_pWgt[i];
      wsum += w*x;
   }
   return wsum;   
} /* end GetInitialVal() */

/******************************************************************************
TiedRespVarWsum::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedRespVarWsum::Write(FILE * pFile, int type)
{
   double val;
   int i;

   val = GetCurrentVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Response Variable (weight)\n");
      for(i = 0; i < m_Num; i++)
      {
         fprintf(pFile, "%s  (%lf)\n", m_pList[i]->GetName(), m_pWgt[i]);
      }
      fprintf(pFile, "Value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end TiedRespVarWsum::WriteToFile() */
