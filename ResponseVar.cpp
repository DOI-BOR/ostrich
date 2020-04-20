/******************************************************************************
File      : ResponseVar.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a single response variable, the optimization analog of regression 
observations.

Version History
05-10-04    lsm   created
01-10-05    lsm   Modified to inherit from abstract response variables (RespVarABC)
******************************************************************************/
#include <stdio.h>
#include <string.h>

#include "ResponseVar.h"
#include "ParameterABC.h"
#include "TiedParamABC.h"

#include "Exception.h"

/******************************************************************************
GetName()

Returns the name of the response variable.
******************************************************************************/
UnchangeableString ResponseVar::GetName(void)
{
  return m_Name;
} /* end ResponseVar() */

/******************************************************************************
CTOR

Dummy constructor of the class.
******************************************************************************/
ResponseVar::ResponseVar(void)
{
   m_Name = NULL;
   m_FileName = NULL;
   m_Keyword = NULL;
   m_Tok = ' ';
   m_bAug = false;
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()
******************************************************************************/
void ResponseVar::Destroy(void)
{
   if(m_Name != NULL)
   {
      delete [] m_Name;
      delete [] m_FileName;
      delete [] m_Keyword;
   } /* end if() */

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
SetCurrentVal()

Sets the current value to the given value.
******************************************************************************/
void ResponseVar::SetCurrentVal(double curVal)
{
  m_CurrentVal = curVal;
} /* end SetCurrentVal() */

/******************************************************************************
SetInitialVal()

Sets the initial value to the given value.
******************************************************************************/
void ResponseVar::SetInitialVal(double initVal)
{
  m_InitialVal = initVal;
} /* end SetInitialVal() */

/******************************************************************************
GetFileName()

Returns the file name associated with the response variable. 
******************************************************************************/
UnchangeableString ResponseVar::GetFileName(void)
{
  return m_FileName;
} /* end GetFileName() */

/******************************************************************************
GetKeyword()

Returns the key word associated with the response variable. The extraction of 
the response variable value depends on the key word as the extracting method 
first looks for the keyword.
******************************************************************************/
UnchangeableString ResponseVar::GetKeyword(void)
{
  return m_Keyword;
} /* end GetKeyword() */

/******************************************************************************
GetLine()

Returns the line number associated with the response variable.
******************************************************************************/
int ResponseVar::GetLine(void)
{
   if(m_pLine != NULL){ m_Line = (int)(m_pLine->GetTransformedVal());}
   if(m_pTiedLine != NULL){ m_Line = (int)(m_pTiedLine->GetEstVal());}
   return m_Line;
} /* end GetLine() */

/******************************************************************************
GetColumn()

Returns the column associated with the response variable
******************************************************************************/
int ResponseVar::GetColumn(void)
{
   if(m_pCol != NULL){ m_Column = (int)(m_pCol->GetTransformedVal());}
   if(m_pTiedCol != NULL){ m_Column= (int)(m_pTiedCol->GetEstVal());}

  return m_Column;
} /* end GetColumn() */

/******************************************************************************
GetCurrentVal()

Returns the current value for the response variable.
******************************************************************************/
double ResponseVar::GetCurrentVal(void)
{
  return m_CurrentVal;
} /* end GetCurrentVal() */

/******************************************************************************
GetInitialVal()

Returns the initial value of the response variable.
******************************************************************************/
double ResponseVar::GetInitialVal(void)
{
  return m_InitialVal;
} /* end GetInitialVal() */

/******************************************************************************
CTOR 

Constructor for the class, using constant values for line and column
******************************************************************************/
ResponseVar::ResponseVar
(
   IroncladString name,
	IroncladString fileName, 
   IroncladString keyword, 
   int line,
	int column,
   char tok,
   bool bAug
)
{   
   int len;

   len = (int)strlen(name) + 1;
   NEW_PRINT("char", len);   
   m_Name = new char[len];   

   len = (int)strlen(keyword) + 1;
   NEW_PRINT("char", len);   
   m_Keyword = new char[len];

   len = (int)strlen(fileName) + 1;
   NEW_PRINT("char", len);   
   m_FileName = new char[len];
   MEM_CHECK(m_FileName);

   strcpy(m_Name, name);
   strcpy(m_Keyword, keyword);
   strcpy(m_FileName, fileName);

   m_InitialVal  = 0.00;
   m_CurrentVal  = 0.00;
   m_Line        = line;
   m_Column      = column;
   m_Tok         = tok;
   m_bAug        = bAug;

   m_pLine       = NULL;
   m_pCol        = NULL;
   m_pTiedLine   = NULL;
   m_pTiedCol    = NULL;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void ResponseVar::Write(FILE * pFile, int type)
{
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         fprintf(pFile, "%E  %E  ", m_InitialVal, m_CurrentVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
	      fprintf(pFile, "%.6lf  %.6lf  ", m_InitialVal, m_CurrentVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_BNR) :
      {
	      fprintf(pFile, "%-12s  initial       current       ", m_Name);
         break;
      }/* end case(WRITE_SCI) */
      default:
      case(WRITE_DBG) :
      {
         m_Line = GetLine();
         m_Column = GetColumn();
         fprintf(pFile, "------Response Variable------\n");
         fprintf(pFile, "Name     : %s\n", m_Name);
         fprintf(pFile, "Filename : %s\n", m_FileName);
         fprintf(pFile, "Keyword  : %s\n", m_Keyword);
         fprintf(pFile, "Line     : %4d    Column : %4d\n", m_Line, m_Column);
         fprintf(pFile, "Token    : %c\n", m_Tok);

         fprintf(pFile, "Initial  : %.6lf  Current : %.6lf\n", m_InitialVal, 
                 m_CurrentVal);
         if(m_pLine != NULL)
         {
            fprintf(pFile, "Line derived from paramter: %s\n", 
            m_pLine->GetName());
         }
         if(m_pCol != NULL)
         {
            fprintf(pFile, "Column derived from paramter: %s\n", 
            m_pCol->GetName());
         }
         if(m_pTiedLine != NULL)
         {
            fprintf(pFile, "Line derived from tied paramter: %s\n", 
            m_pTiedLine->GetName());
         }
         if(m_pTiedCol != NULL)
         {
            fprintf(pFile, "Column derived from tied paramter: %s\n", 
            m_pTiedCol->GetName());
         }
         break;
      }/* end case(WRITE_SCI) */
   }/* end switch() */
} /* end Write() */

/******************************************************************************
WriteSim()

Writes simulated output to pFile.
******************************************************************************/
void ResponseVar::WriteSim(FILE * pFile, int type)
{
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         fprintf(pFile, "%E  ", m_CurrentVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
	      fprintf(pFile, "%.6lf  ", m_CurrentVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_BNR) :
      {
	      fprintf(pFile, "%-12s  ", m_Name);
         break;
      }/* end case(WRITE_SCI) */
      default:
      case(WRITE_DBG) :
      {
         fprintf(pFile, "%-12s = %.6lf\n", m_Name, m_CurrentVal);
         break;
      }/* end case(WRITE_DBG) */
   }/* end switch() */
} /* end WriteSim() */
