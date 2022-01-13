/******************************************************************************
File      : TiedParam.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Definition of various 'tied' parameters. Tied parameters are variables in the 
model which are computed from the values of one or more model parameters. The 
ABC for the tied parameters encapsulates the interface used by other Ostrich 
modules, allowing various specific tied parameter relationships (linear, 
exponential, etc.) to be implemented as needed with minimal code change (just 
need to add the specific tied parameter class and some additional input file
parsing).

These specific tied-parameter classes are supported:
TiedParamLin1         : linear function of one parameter
TiedParamLin2         : linear function of two parameters
TiedParamExp          : exponential function of one parameter
TiedParamLog          : logarithmic function of one parameter
TiedDistXY            : distance between two (x,y) parameters
TiedSimpleParamRatio  : ratio of two parameters (ax + b) / (cy + d)
TiedComplexParamRatio : complex ratio of three parameters
      (a*x*y*z + b*x*y + c*x*z + d*y*z + e*x + f*y + g*z + h) / 
      (i*x*y*z + j*x*y + k*x*z + l*y*z + m*x + n*y + o*z + p)
TiedParamConstant     : parameter is assigned a constant value
TiedParamWsum  : weighted sum of one or more parameters

Version History
07-07-04    lsm   Created
09-13-04    lsm   added distance (TiedDistXY)
01-01-07    lsm   added ratios (TiedParamSimpleRatio and TiedParamComplexRatio)
08-06-09    lsm   added constants
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "TiedParamABC.h"
#include "ParameterABC.h"

#include "FortranSupportUtilities.h"
#include "Exception.h"
#include "Utility.h"

UnchangeableString GetMetaName(MetaParameter * mp)
{
   if(mp == NULL) return NULL;
   if(mp->pParam == NULL) return NULL;
   if(mp->type == RGLR_PARAMETER) return (((ParameterABC *)(mp->pParam))->GetName());
   return (((TiedParamABC *)(mp->pParam))->GetName());
   return NULL;
}/* end GetMetaName() */

double GetMetaVal(MetaParameter * mp)
{
   if(mp == NULL) return 0.00;
   if(mp->pParam == NULL) return 0.00;
   if(mp->type == RGLR_PARAMETER) return (((ParameterABC *)(mp->pParam))->GetTransformedVal());
   return (((TiedParamABC *)(mp->pParam))->GetEstVal());
   return 0.00;
}/* end GetMetaVal() */

/******************************************************************************
TiedParamLin1::Destroy()
******************************************************************************/
void TiedParamLin1::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamLin1)
******************************************************************************/
TiedParamLin1::TiedParamLin1(void)
{
   m_pName = NULL;
   m_Tie.pParam = NULL;
   m_pFixFmt = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamLin1)
******************************************************************************/
TiedParamLin1::TiedParamLin1(IroncladString name, MetaParameter * p1, 
                             UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   m_pFixFmt = new char[len];
   MEM_CHECK(m_pFixFmt);

   strcpy(m_pName, name);

   m_Tie.pParam = p1->pParam;  
   m_Tie.type = p1->type;

   //parse the config string to determine value for C0 and C1
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin1()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);
   //extract fortran formatting
   pTok += j;
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   char valStr[DEF_STR_SZ];
   bool bOk;
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamLin1::GetEstVal(void)
{
   double y, x;
   x = GetMetaVal(&m_Tie);
   y = m_C1*x + m_C0;
   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamLin1::GetValAsStr()
******************************************************************************/
void TiedParamLin1::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamLin1::GetValAsStr() */

/******************************************************************************
TiedParamLin1::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamLin1::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied Param = %s\n", GetMetaName(&m_Tie));
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
} /* end TiedParamLin1::Write() */

/******************************************************************************
TiedParamLin2::Destroy()
******************************************************************************/
void TiedParamLin2::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamLin2)
******************************************************************************/
TiedParamLin2::TiedParamLin2(void)
{
   m_pName = NULL;
   m_pFixFmt = NULL;
   m_Tie1.pParam = NULL;
   m_Tie2.pParam = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamLin2)
******************************************************************************/
TiedParamLin2::TiedParamLin2(IroncladString name, MetaParameter * p1, 
                             MetaParameter * p2, UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];
   char valStr[DEF_STR_SZ];
   bool bOk;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   m_pFixFmt = new char[len];
   MEM_CHECK(m_pFixFmt);

   strcpy(m_pName, name);

   m_Tie1.pParam = p1->pParam;
   m_Tie2.pParam = p2->pParam;
   m_Tie1.type = p1->type;
   m_Tie2.type = p2->type;

   //parse the config string to determin valuese for C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin2()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin2()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin2()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);
   //extract fortran formatting
   pTok += j;
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedParamLin2::GetValAsStr()
******************************************************************************/
void TiedParamLin2::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamLin2::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamLin2::GetEstVal(void)
{
   double y, x1, x2;
   x1 = GetMetaVal(&m_Tie1);
   x2 = GetMetaVal(&m_Tie2);
   y = m_C3*x1*x2 + m_C2*x2 + m_C1*x1 + m_C0;
   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamLin2::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamLin2::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied Param #1 = %s\n", GetMetaName(&m_Tie1));
      fprintf(pFile, "Tied Param #2 = %s\n", GetMetaName(&m_Tie2));
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
} /* end TiedParamLin2::Write() */

/******************************************************************************
TiedParamExp::Destroy()
******************************************************************************/
void TiedParamExp::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;   
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamExp)
******************************************************************************/
TiedParamExp::TiedParamExp(void)
{
   m_pName = NULL;
   m_Tie.pParam = NULL;
   m_pFixFmt = NULL;
   m_Base = 0.00;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamExp)
******************************************************************************/
TiedParamExp::TiedParamExp(IroncladString name, MetaParameter * p1, 
                           UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];
   char valStr[DEF_STR_SZ];
   bool bOk;
   
   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   m_pFixFmt = new char[len];
   MEM_CHECK(m_pFixFmt);

   strcpy(m_pName, name);

   m_Tie.pParam = p1->pParam;
   m_Tie.type = p1->type;

   //parse the config string to determine values for base, C0, C1, and C2
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamExp()");
   if(strcmp(tmpStr, "exp") == 0){ m_Base = 2.718;}
   else{ m_Base = atof(tmpStr);}
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamExp()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamExp()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);
   //extract fortran formatting
   pTok += j;
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamExp::GetEstVal(void)
{
   double y, x;

   x = GetMetaVal(&m_Tie);

   if(m_Base == 2.718) //use exp
   {
      y = (m_C2 * exp(m_C1*x)) + m_C0;
   }
   else //use pow
   {
      y = (m_C2 * pow(m_Base, (m_C1*x))) + m_C0;
   }   

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamExp::GetValAsStr()
******************************************************************************/
void TiedParamExp::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamExp::GetValAsStr() */

/******************************************************************************
TiedParamExp::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamExp::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied Param = %s\n", GetMetaName(&m_Tie));
      fprintf(pFile, "Exponent Base = %lf\n", m_Base);
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
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
} /* end TiedParamExp::Write() */

/******************************************************************************
TiedParamLog::Destroy()
******************************************************************************/
void TiedParamLog::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamLog)
******************************************************************************/
TiedParamLog::TiedParamLog(void)
{
   m_pName = NULL;
   m_Tie.pParam = NULL;
   m_pFixFmt = NULL;
   m_Base = 0.00;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamLog)
******************************************************************************/
TiedParamLog::TiedParamLog(IroncladString name, MetaParameter * p1, 
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

   m_Tie.pParam = p1->pParam;
   m_Tie.type = p1->type;

   //parse the config string to determine values for base, C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   if(strcmp(tmpStr, "ln") == 0){ m_Base = 2.718;}
   else{ m_Base = atof(tmpStr);}
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);
   //extract fortran formatting
   m_pFixFmt = new char[len];
   pTok += j;   
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   char valStr[DEF_STR_SZ];
   bool bOk;
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedParamLog::GetValAsStr()
******************************************************************************/
void TiedParamLog::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamLog::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamLog::GetEstVal(void)
{
   double y, x, N;

   x = GetMetaVal(&m_Tie);
   N = (m_C2*x) + m_C1;

   if(m_Base == 2.718) //use exp
   {
      y = (m_C3 * log(N)) + m_C0;
   }
   else if(m_Base == 10.00) //use log10
   {
      y = (m_C3 * log10(N)) + m_C0;
   }   
   else //convert from log10, using loga(N) = log10(N)/log10(a)
   {
      y = (m_C3 * (log10(N)/log10(m_Base))) + m_C0;
   }

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamLog::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamLog::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied Param = %s\n", GetMetaName(&m_Tie));
      fprintf(pFile, "Log Base = %lf\n", m_Base);
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
} /* end TiedParamLog::Write() */

/******************************************************************************
TiedDistXY::Destroy()
******************************************************************************/
void TiedDistXY::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedDistXY)
******************************************************************************/
TiedDistXY::TiedDistXY(void)
{
   m_pName = NULL;
   m_pFixFmt = NULL;
   m_X1.pParam = NULL;
   m_Y1.pParam = NULL;
   m_X2.pParam = NULL;
   m_Y2.pParam = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedDistXY)
******************************************************************************/
TiedDistXY::TiedDistXY(IroncladString name, MetaParameter * px1, 
                       MetaParameter * py1,  MetaParameter * px2, 
                       MetaParameter * py2, UnmoveableString configStr)
{
   int len;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_X1.pParam = px1->pParam;
   m_Y1.pParam = py1->pParam;
   m_X2.pParam = px2->pParam;
   m_Y2.pParam = py2->pParam;
   m_X1.type = px1->type;
   m_Y1.type = py1->type;
   m_X2.type = px2->type;
   m_Y2.type = py2->type;

   //extract fortran formatting
   m_pFixFmt = new char[len];
   char * pTok = configStr;
   if(pTok != NULL)
     strcpy(m_pFixFmt, pTok);
   else
     strcpy(m_pFixFmt, "free");
   MyTrim(m_pFixFmt);
   //check fixed format setting
   char valStr[DEF_STR_SZ];
   bool bOk;
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedDistXY::GetValAsStr()
******************************************************************************/
void TiedDistXY::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedDistXY::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedDistXY::GetEstVal(void)
{
   double d, x1, x2, y1, y2;

   x1 = GetMetaVal(&m_X1);
   x2 = GetMetaVal(&m_X2);
   y1 = GetMetaVal(&m_Y1);
   y2 = GetMetaVal(&m_Y2);

   d = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

   return d;
} /* end GetEstVal() */

/******************************************************************************
TiedDistXY::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedDistXY::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied X1 = %s\n", GetMetaName(&m_X1));
      fprintf(pFile, "Tied Y1 = %s\n", GetMetaName(&m_Y1));
      fprintf(pFile, "Tied X2 = %s\n", GetMetaName(&m_X2));
      fprintf(pFile, "Tied Y2 = %s\n", GetMetaName(&m_Y2));
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
} /* end TiedDistXY::Write() */

/******************************************************************************
TiedParamSimpleRatio::Destroy()
******************************************************************************/
void TiedParamSimpleRatio::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamSimpleRatio)
******************************************************************************/
TiedParamSimpleRatio::TiedParamSimpleRatio(void)
{
   m_pName = NULL;
   m_Tie1.pParam = NULL;
   m_Tie2.pParam = NULL;
   m_pFixFmt = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamSimpleRatio)
******************************************************************************/
TiedParamSimpleRatio::TiedParamSimpleRatio(IroncladString name, MetaParameter * p1, 
                               MetaParameter * p2, UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_Tie1.pParam = p1->pParam;
   m_Tie2.pParam = p2->pParam;
   m_Tie1.type = p1->type;
   m_Tie2.type = p2->type;

   //parse the config string to determine values for C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamSimpleRatio()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamSimpleRatio()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamSimpleRatio()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);
   //extract fortran formatting
   m_pFixFmt = new char[len];
   pTok += j;   
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   char valStr[DEF_STR_SZ];
   bool bOk;
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedParamSimpleRatio::GetValAsStr()
******************************************************************************/
void TiedParamSimpleRatio::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamSimpleRatio::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamSimpleRatio::GetEstVal(void)
{
   double y, x1, x2;

   x1 = GetMetaVal(&m_Tie1);
   x2 = GetMetaVal(&m_Tie2);
   y = (m_C3*x1 + m_C2) / (m_C1*x2 + m_C0);

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamSimpleRatio::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamSimpleRatio::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied Param #1 = %s\n", GetMetaName(&m_Tie1));
      fprintf(pFile, "Tied Param #2 = %s\n", GetMetaName(&m_Tie2));
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
      fprintf(pFile, "C3 = %lf\n", m_C3);
      fprintf(pFile, "value = %lf\n", val);
      fprintf(pFile, "function = (C3*P2 + C2)/(C1*P1 + C0)\n");
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
} /* end TiedParamSimpleRatio::Write() */

/******************************************************************************
TiedParamComplexRatio::Destroy()
******************************************************************************/
void TiedParamComplexRatio::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamComplexRatio)
******************************************************************************/
TiedParamComplexRatio::TiedParamComplexRatio(void)
{
   m_pName = NULL;
   m_pFixFmt = NULL;
   m_X.pParam = NULL;
   m_Y.pParam = NULL;
   m_Z.pParam = NULL;
   for(int i = 0; i < 8; i++){ m_N[i] = m_D[i] = 0.00;}
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamComplexRatio)
******************************************************************************/
TiedParamComplexRatio::TiedParamComplexRatio
(
   IroncladString name, 
   MetaParameter * p1, 
   MetaParameter * p2, 
   MetaParameter * p3, 
   UnmoveableString configStr
)
{
   int i, j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_X.pParam = p1->pParam;
   m_Y.pParam = p2->pParam;
   m_Z.pParam = p3->pParam;
   m_X.type = p1->type;
   m_Y.type = p2->type;
   m_Z.type = p3->type;

   //parse the config string to determine values for numer and denom coeffs
   pTok = configStr;
   for(i = 0; i < 8; i++)
   {
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, i, 8, "TiedParamSimpleRatio()");
      m_N[7-i] = atof(tmpStr);
      pTok += j;
   }
   for(i = 0; i < 8; i++)
   {
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, i, 8, "TiedParamSimpleRatio()");
      m_D[7-i] = atof(tmpStr);
      pTok += j;
   }
   //extract fortran formatting
   m_pFixFmt = new char[len];
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   char valStr[DEF_STR_SZ];
   bool bOk;
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedParamComplexRatio::GetValAsStr()
******************************************************************************/
void TiedParamComplexRatio::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamComplexRatio::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamComplexRatio::GetEstVal(void)
{
   double x, y, z, num, den, val;

   x = GetMetaVal(&m_X);
   y = GetMetaVal(&m_Y);
   z = GetMetaVal(&m_Z);

   num = m_N[7]*x*y*z + m_N[6]*x*y + m_N[5]*x*z + m_N[4]*y*z + 
         m_N[3]*x     + m_N[2]*y   + m_N[1]*z   + m_N[0];

   den = m_D[7]*x*y*z + m_D[6]*x*y + m_D[5]*x*z + m_D[4]*y*z + 
         m_D[3]*x     + m_D[2]*y   + m_D[1]*z   + m_D[0];

   val = (num / den);
   return val;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamComplexRatio::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamComplexRatio::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

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
      fprintf(pFile, "Tied Param #1 (X) = %s\n", GetMetaName(&m_X));
      fprintf(pFile, "Tied Param #2 (Y) = %s\n", GetMetaName(&m_Y));
      fprintf(pFile, "Tied Param #3 (Z) = %s\n", GetMetaName(&m_Z));
      fprintf(pFile, "A = %lf\n", m_N[7]);
      fprintf(pFile, "B = %lf\n", m_N[6]);
      fprintf(pFile, "C = %lf\n", m_N[5]);
      fprintf(pFile, "D = %lf\n", m_N[4]);
      fprintf(pFile, "E = %lf\n", m_N[3]);
      fprintf(pFile, "F = %lf\n", m_N[2]);
      fprintf(pFile, "G = %lf\n", m_N[1]);
      fprintf(pFile, "H = %lf\n", m_N[0]);
      fprintf(pFile, "I = %lf\n", m_D[7]);
      fprintf(pFile, "J = %lf\n", m_D[6]);
      fprintf(pFile, "K = %lf\n", m_D[5]);
      fprintf(pFile, "L = %lf\n", m_D[4]);
      fprintf(pFile, "M = %lf\n", m_D[3]);
      fprintf(pFile, "N = %lf\n", m_D[2]);
      fprintf(pFile, "O = %lf\n", m_D[1]);
      fprintf(pFile, "P = %lf\n", m_D[0]);
      fprintf(pFile, "value = %lf\n", val);
      fprintf(pFile, "function = (Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H) /");
      fprintf(pFile, "           (Ixyz + Jxy + Kxz + Lyz + Mx + Ny + Oz + P)\n");
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
} /* end TiedParamComplexRatio::Write() */

/******************************************************************************
TiedParamConstant::Destroy()
******************************************************************************/
void TiedParamConstant::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   delete [] m_pStringLiteral;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamConstant)
******************************************************************************/
TiedParamConstant::TiedParamConstant(void)
{
   m_pName = NULL;
   m_pFixFmt = NULL;
   m_pStringLiteral = NULL;
   m_val = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamConstant)
******************************************************************************/
TiedParamConstant::TiedParamConstant(IroncladString name, UnmoveableString pVal)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   //parse the config string to determine value of constant
   pTok = pVal;
   j = ExtractString(pTok, tmpStr);
   //check for string literal
   if(((tmpStr[0] == '"') && (tmpStr[strlen(tmpStr)-1] == '"')) ||
      ((tmpStr[0] == '\'') && (tmpStr[strlen(tmpStr) - 1] == '\'')))
   {
      m_pStringLiteral = new char[strlen(tmpStr)];
      strcpy(m_pStringLiteral, &(tmpStr[1]));
      m_pStringLiteral[strlen(tmpStr) - 2] = NULLSTR;
   }
   else
   {
      m_val = atof(tmpStr);
      //extract fortran formatting
      m_pFixFmt = new char[len];
      pTok += j;   
      strcpy(m_pFixFmt, pTok);
      MyTrim(m_pFixFmt);
      //check fixed format setting
      char valStr[DEF_STR_SZ];
      bool bOk;
      bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
      if(!bOk) strcpy(m_pFixFmt, "free");
   }
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedParamConstant::GetValAsStr()
******************************************************************************/
void TiedParamConstant::GetValAsStr(UnmoveableString valStr)
{
  if(m_pStringLiteral != NULL)
  {
     strcpy(valStr, m_pStringLiteral);
  }
  else
  {
      bool bOk;
      double val = GetEstVal();
      char * fmt = m_pFixFmt;
      if(strcmp(fmt, "free") == 0)
      {
         GetPreciseValAsStr(valStr, val);
      }
      else
      {
         bOk = GetFixedFormatValAsStr(valStr, val, fmt);
         if(!bOk)
         {
            GetPreciseValAsStr(valStr, val);
         }
      }
   }
}/* end TiedParamConstant::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamConstant::GetEstVal(void)
{
   return m_val;
} /* end GetEstVal() */

/******************************************************************************
TiedParamConstant::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamConstant::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      if(m_pStringLiteral != NULL)
      {
         fprintf(pFile, "%-12s  ", m_pStringLiteral);
      }
      else
      {
         fprintf(pFile,"%E  ", val);
      }
   }
   else if (type == WRITE_DEC)
   {
      if (m_pStringLiteral != NULL)
      {
         fprintf(pFile, "%-12s  ", m_pStringLiteral);
      }
      else
      {
         fprintf(pFile,"%.6lf  ", val);
      }
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      if (m_pStringLiteral != NULL)
      {
         fprintf(pFile, "Value = %s\n", m_pStringLiteral);
      }
      else
      {
         fprintf(pFile, "Value = %lf\n", val);
      }
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      if (m_pStringLiteral != NULL)
      {
         fprintf(pFile, "%-18s : %s\n", m_pName, m_pStringLiteral);
      }
      else
      {
         fprintf(pFile, "%-18s : %E\n", m_pName, val);
      }
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end TiedParamConstant::Write() */

/******************************************************************************
TiedParamWsum::Destroy()
******************************************************************************/
void TiedParamWsum::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pTie;
   delete [] m_pWgt;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamWsum)
******************************************************************************/
TiedParamWsum::TiedParamWsum(void)
{
   m_pName = NULL;
   m_pWgt  = NULL;
   m_pTie = NULL;
   m_pFixFmt = NULL;
   m_Num = 0;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamWsum)
******************************************************************************/
TiedParamWsum::TiedParamWsum(IroncladString name, MetaParameter * pList, int num,
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

   NEW_PRINT("MetaParameter", num);
   m_pTie = new MetaParameter[num];
   MEM_CHECK(m_pTie);
   for(i = 0; i < num; i++)
   {
      m_pTie[i].pParam = pList[i].pParam;
      m_pTie[i].type = pList[i].type;
   }

   m_Num = num;

   NEW_PRINT("double", num);
   m_pWgt = new double[num];
   MEM_CHECK(m_pWgt);
   for(i = 0; i < num; i++) m_pWgt[i] = 0.00;

   //parse the config string to determine value for weights
   pTok = configStr;
   for(i = 0; i < num; i++)
   {
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, i, num, "TiedParamWsum()");
      m_pWgt[i] = atof(tmpStr);
      pTok += j;
   }

   //extract fortran formatting
   m_pFixFmt = new char[len];
   strcpy(m_pFixFmt, pTok);
   MyTrim(m_pFixFmt);
   //check fixed format setting
   char valStr[DEF_STR_SZ];
   bool bOk;
   bOk = GetFixedFormatValAsStr(valStr, 0.00, m_pFixFmt);
   if(!bOk) strcpy(m_pFixFmt, "free");
      
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
TiedParamWsum::GetValAsStr()
******************************************************************************/
void TiedParamWsum::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetEstVal();
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, val);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, val, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, val);
      }
   }
}/* end TiedParamWsum::GetValAsStr() */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamWsum::GetEstVal(void)
{
   int i;
   double wsum, val, wgt;

   wsum = 0.00;

   for(i = 0; i < m_Num; i++)
   {
      val = GetMetaVal(&(m_pTie[i]));
      wgt = m_pWgt[i];
      wsum += (val*wgt);
   }
   return wsum;
} /* end GetEstVal() */

/******************************************************************************
TiedParamWsum::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamWsum::Write(FILE * pFile, int type)
{
   int i;
   double val;

   val = GetEstVal();

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
      for(i = 0; i < m_Num; i++)
      {
         fprintf(pFile, "Tied Param #%d = %s\n", (i+1), GetMetaName(&(m_pTie[i])));
         fprintf(pFile, "Weight #%d = %E\n", (i+1), m_pWgt[i]);
      }
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
} /* end TiedParamWsum::Write() */
