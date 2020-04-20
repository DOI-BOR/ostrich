/******************************************************************************
File      : ComboParam.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

A series of classes that encapsulate various types of combinatorial parameters.

Version History
03-18-04    lsm   created
07-08-04    lsm   integrated with ParameterABC class, added WRITE_OPT support
08-17-04    lsm   RAM fragmentation fixes, inherits from ParameterABC
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ParameterABC.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
ComboIntParam::Destroy()
******************************************************************************/
void ComboIntParam::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pCombos;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR
******************************************************************************/
ComboIntParam::ComboIntParam(void)
{
   m_pName = NULL;
   m_pCombos = NULL;
   m_CurIdx = -1;
   m_NumCombos = 0;
   m_InitIdx = -1;

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR
******************************************************************************/
ComboIntParam::ComboIntParam(IroncladString name, UnmoveableString configStr)
{
   int i, j, init, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];
   m_CurIdx = -1;
   m_NumCombos = 0;
   m_InitIdx = -1;

   len = (int)strlen(name) + 1;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);
   strcpy(m_pName, name);

   pTok = configStr;
   while(*pTok == '\t'){ pTok++;}
   //extract initial value
   j = ExtractColString(pTok, tmpStr, '\t');
   j = ValidateExtraction(j, 1, 1, "ComboIntParam()");
   init = atoi(tmpStr);
   pTok += j;
   //extract number of combinations
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "ComboIntParam()");

   m_NumCombos = atoi(tmpStr);
   if(m_NumCombos <= 0)
   {
      LogError(ERR_FILE_IO, "ComboIntParam(): Invalid number of combinations");
      ExitProgram(1);
   }
   pTok += j;
   NEW_PRINT("int", m_NumCombos);
   m_pCombos = new int[m_NumCombos];   
   MEM_CHECK(m_pCombos);   
   //extract the combinations
   while(*pTok == '\t'){ pTok++;}
   for(i = 0; i < m_NumCombos; i++)
   {
      j = ExtractColString(pTok, tmpStr, '\t');
      j = ValidateExtraction(j, i, m_NumCombos, "ComboIntParam()");
      m_pCombos[i] = atoi(tmpStr);
      pTok += j;
      if(m_pCombos[i] == init){ m_InitIdx = m_CurIdx = i;}
   }

   if((m_InitIdx == -1) || (m_CurIdx == -1))
   {
      LogError(ERR_FILE_IO, "ComboIntParam(): Invalid initial parameter value");
      ExitProgram(1);
   }/* end if() */
      
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SetEstVal()

Set the estimated value of the parameter.
******************************************************************************/
double ComboIntParam::SetEstVal(double Idx)
{ 
   int tmp = (int)(Idx + 0.5); //round up real value

   if((tmp < m_NumCombos) && (tmp >= 0)){ m_CurIdx = tmp;}

   return 0.00;
}/* end SetEstVal() */

/******************************************************************************
Write()

Writes formatted output to pFile argument.
******************************************************************************/
void ComboIntParam::Write(FILE * pFile, int type)
{
   int val, i;

   val = m_pCombos[m_CurIdx];

   if((type == WRITE_SCI) || (type == WRITE_DEC))
   {
      fprintf(pFile,"%-12d  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Initial Value   (%d) %d\n", m_InitIdx, m_pCombos[m_InitIdx]);
      fprintf(pFile, "Estimated Value (%d) %d\n", m_CurIdx, val);
      fprintf(pFile, "Lower Bound 0\n");
      fprintf(pFile, "Upper Bound %d\n", (m_NumCombos-1));
      fprintf(pFile, "Possible Values\n");
      for(i = 0; i < m_NumCombos; i++)
      {
         fprintf(pFile, "(%d) %d\n", i, m_pCombos[i]);
      }/* end for() */
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %d\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end writeToFile() */

/******************************************************************************
GetValAsStr()
******************************************************************************/
void ComboDblParam::GetValAsStr(UnmoveableString valStr)
{
	GetPreciseValAsStr(valStr, m_pCombos[m_CurIdx]); 
}/* end GetValAsStr() */

/******************************************************************************
Destroy()
******************************************************************************/
void ComboDblParam::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pCombos;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR
******************************************************************************/
ComboDblParam::ComboDblParam(void)
{
   m_pName = NULL;
   m_pCombos = NULL;
   m_CurIdx = -1;
   m_NumCombos = 0;
   m_InitIdx = -1;

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR
******************************************************************************/
ComboDblParam::ComboDblParam(IroncladString name, UnmoveableString configStr)
{
   int i, j, len;
   double init;
   char * pTok;
   char tmpStr[DEF_STR_SZ];
   m_CurIdx = -1;
   m_NumCombos = 0;
   m_InitIdx = -1;

   len = (int)strlen(name) + 1;
   NEW_PRINT("char", len);
   m_pName = new char[len];   
   MEM_CHECK(m_pName);
   strcpy(m_pName, name);

   pTok = configStr;
   while(*pTok == '\t'){ pTok++;}
   //extract initial value
   j = ExtractColString(pTok, tmpStr, '\t');
   j = ValidateExtraction(j, 1, 1, "ComboDblParam()");
   init = atof(tmpStr);
   pTok += j;
   //extract number of cominations
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "ComboDblParam()");
   m_NumCombos = atoi(tmpStr);
   if(m_NumCombos <= 0)
   {
      LogError(ERR_FILE_IO, "ComboDblParam(): Invalid number of combinations");
      ExitProgram(1);
   }
   pTok += j;
   NEW_PRINT("double", m_NumCombos);
   m_pCombos = new double[m_NumCombos];   
   MEM_CHECK(m_pCombos);
   //extract the combinations
   while(*pTok == '\t'){ pTok++;}
   for(i = 0; i < m_NumCombos; i++)
   {
      j = ExtractColString(pTok, tmpStr, '\t');
      j = ValidateExtraction(j, i, m_NumCombos, "ComboDblParam()");
      m_pCombos[i] = atof(tmpStr);
      pTok += j;
      if(m_pCombos[i] == init){ m_InitIdx = m_CurIdx = i;}
   }

   if((m_InitIdx == -1) || (m_CurIdx == -1))
   {
      LogError(ERR_FILE_IO, "ComboDblParam(): Invalid initial parameter value");
      ExitProgram(1);
   }/* end if() */
               
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SetEstVal()

Set the estimated value of the parameter.
******************************************************************************/
double ComboDblParam::SetEstVal(double Idx)
{ 
   int tmp = (int)(Idx + 0.5); //round up real value

   if((tmp < m_NumCombos) && (tmp >= 0)){ m_CurIdx = tmp;}

   return 0.00;
}/* end SetEstVal() */

/******************************************************************************
Write()

Writes formatted output to pFile argument.
******************************************************************************/
void ComboDblParam::Write(FILE * pFile, int type)
{
   double val;
   int i;

   val = m_pCombos[m_CurIdx];

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
      fprintf(pFile, "Name %s\n", m_pName);
      fprintf(pFile, "Initial Value   (%d) %E\n", m_InitIdx, m_pCombos[m_InitIdx]);
      fprintf(pFile, "Estimated Value (%d) %E\n", m_CurIdx, val);
      fprintf(pFile, "Lower Bound 0\n");
      fprintf(pFile, "Upper Bound %d\n", (m_NumCombos-1));
      fprintf(pFile, "Possible Values\n");
      for(i = 0; i < m_NumCombos; i++)
      {
         fprintf(pFile, "(%d) %E\n", i, m_pCombos[i]);
      }/* end for() */
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
} /* end writeToFile() */

/******************************************************************************
Destroy()
******************************************************************************/
void ComboStrParam::Destroy(void)
{
   delete [] m_pName;
   
   for(int i = 0; i < m_NumCombos; i++){delete [] m_pCombos[i];}
   delete [] m_pCombos;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR
******************************************************************************/
ComboStrParam::ComboStrParam(void)
{
   m_pName = NULL;
   m_pCombos = NULL;
   m_CurIdx = -1;
   m_NumCombos = 0;
   m_InitIdx = -1;

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR
******************************************************************************/
ComboStrParam::ComboStrParam(IroncladString name, UnmoveableString configStr)
{
   int i, j, len;
   char * init;
   char * pTok;
   char tmpStr[DEF_STR_SZ];
   m_CurIdx = -1;
   m_NumCombos = 0;
   m_InitIdx = -1;

   len = (int)strlen(name) + 1;
   NEW_PRINT("char", len);
   m_pName = new char[len];   
   MEM_CHECK(m_pName);
   strcpy(m_pName, name);

   pTok = configStr;
   while(*pTok == '\t'){ pTok++;}
   //extract initial value
   j = ExtractColString(pTok, tmpStr, '\t');
   j = ValidateExtraction(j, 1, 1, "ComboStrParam()");
   len = (int)strlen(tmpStr) + 1;
   NEW_PRINT("char", len);
   init = new char[len];   
   MEM_CHECK(init);
   strcpy(init, tmpStr);
   pTok += j;
   //extract number of cominations
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "ComboStrParam()");
   m_NumCombos = atoi(tmpStr);
   if(m_NumCombos <= 0)
   {
      LogError(ERR_FILE_IO, "ComboStrParam(): Invalid number of combinations");
      ExitProgram(1);
   }
   pTok += j;
   NEW_PRINT("char *", m_NumCombos);
   m_pCombos = new char *[m_NumCombos];   
   MEM_CHECK(m_pCombos);
   //extract the combinations
   while(*pTok == '\t'){ pTok++;}
   for(i = 0; i < m_NumCombos; i++)
   {
      j = ExtractColString(pTok, tmpStr, '\t');
      j = ValidateExtraction(j, i, m_NumCombos, "ComboStrParam()");
      len = (int)strlen(tmpStr) + 1;
      NEW_PRINT("char", len);
      m_pCombos[i] = new char[len];      
      MEM_CHECK(m_pCombos[i]);
      strcpy(m_pCombos[i], tmpStr);
      pTok += j;
      if(strcmp(m_pCombos[i],init) == 0){ m_InitIdx = m_CurIdx = i;}
   }

   if((m_InitIdx == -1) || (m_CurIdx == -1))
   {
      LogError(ERR_FILE_IO, "ComboStrParam(): Invalid initial parameter value");
      ExitProgram(1);
   }/* end if() */

   delete [] init;               
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SetEstVal()

Set the estimated value of the parameter.
******************************************************************************/
double ComboStrParam::SetEstVal(double Idx)
{ 
   int tmp = (int)(Idx + 0.5); //round up real value

   if((tmp < m_NumCombos) && (tmp >= 0)){ m_CurIdx = tmp;}

   return 0.00;
}/* end SetEstVal() */

/******************************************************************************
GetValAsStr()

Get the value of the combinatorial parameter, in a format that can be written
to a template file. This is where OST_NULL gets converted into a null string.
******************************************************************************/
void ComboStrParam::GetValAsStr(UnmoveableString valStr)
{
   if(strcmp(m_pCombos[m_CurIdx], "OST_NULL") == 0){ strcpy(valStr, "");}
   else{ strcpy(valStr, m_pCombos[m_CurIdx]);}
}/* end GetValAsStr() */

/******************************************************************************
Write()

Writes formatted output to pFile argument.
******************************************************************************/
void ComboStrParam::Write(FILE * pFile, int type)
{
   char * val;
   int i;

   val = m_pCombos[m_CurIdx];

   if((type == WRITE_SCI) || (type == WRITE_DEC))
   {
      fprintf(pFile,"%-12s  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name %s\n", m_pName);
      fprintf(pFile, "Initial Value   (%d) %s\n", m_InitIdx, m_pCombos[m_InitIdx]);
      fprintf(pFile, "Estimated Value (%d) %s\n", m_CurIdx, val);      
      fprintf(pFile, "Lower Bound 0\n");
      fprintf(pFile, "Upper Bound %d\n", (m_NumCombos-1));
      fprintf(pFile, "Possible Values\n");
      for(i = 0; i < m_NumCombos; i++)
      {
         fprintf(pFile, "(%d) %s\n", i, m_pCombos[i]);
      }/* end for() */
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %s\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end writeToFile() */
