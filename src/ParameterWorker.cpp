/******************************************************************************
File      : Parameter.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates various types of parameters. Parameters are variables in the model 
which are to be calibrated or optimized. Five parameter classes are defined:
   RealParam : continously varying parameters
   IntParam  : discrete (integer) parameters
   ComboIntParam : integer combinations
   ComboDblParam : real (double) combinations
   ComboStrParam : text (string) combinations

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-25-03    lsm   Modified to support multiple stages of unit conversion.
03-05-04    lsm   replaced magic # with NEARLY_ZERO
07-05-04    lsm   Added ParameterABC to wrap functionality of all param. types
01-05-05    lsm   Replaced truncation with round-up in IntParam::SetEstVal()
03-03-05    lsm   Added support for ON/OFF parameter threshold
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ParameterWorker.h"

#include "Exception.h"
#include "FortranSupportUtilities.h"
#include "Utility.h"

/******************************************************************************
RealParam::Destroy() 
******************************************************************************/
void RealParamWorker::Destroy(void) {
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (RealParam)
******************************************************************************/
RealParamWorker::RealParamWorker(void) {
   m_pName = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (RealParam)
******************************************************************************/
RealParamWorker::RealParamWorker(IroncladString name, double value, IroncladString fixFmt) {
   int len;

   // Set the name into the parameter
   len = (int)strlen(name) + 10;  
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);
   strcpy(m_pName, name);

   // Set the format into the parameter
   len = (int)strlen(fixFmt) + 10;  
   NEW_PRINT("char", len);
   m_pFixFmt = new char[len];
   MEM_CHECK(m_pFixFmt);
   strcpy(m_pFixFmt, fixFmt);

   // Set the value into the parameter
   m_value = value;
      
   IncCtorCount();
} 


/******************************************************************************
RealParam::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void RealParamWorker::Write(FILE * pFile, int type) {

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", m_value);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", m_value);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s  ", m_pName);
      fprintf(pFile, "Value %E\n", m_value);

   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, m_value);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} 

/******************************************************************************
RealParam::GetValAsStr()
******************************************************************************/
void RealParamWorker::GetValAsStr(UnmoveableString valStr) {
   bool bOk;
   char * fmt = m_pFixFmt;
   if(strcmp(fmt, "free") == 0)
   {
      GetPreciseValAsStr(valStr, m_value);
   }
   else
   {
      bOk = GetFixedFormatValAsStr(valStr, m_value, fmt);
      if(!bOk)
      {
         GetPreciseValAsStr(valStr, m_value);
      }
   }
}/* end RealParam::GetValAsStr() */

/******************************************************************************
IntParam::Destroy()
******************************************************************************/
void IntParamWorker::Destroy(void) {
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (IntParam)
******************************************************************************/
IntParamWorker::IntParamWorker(void) {
   m_pName = NULL;
   IncCtorCount();
} /* end default CTOR */



/******************************************************************************
CTOR (IntParam)
******************************************************************************/
IntParamWorker::IntParamWorker(IroncladString name, int initialValue) {
   int len;

   // Set the name into the parameter
   len = (int)strlen(name) + 10;  
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);
   strcpy(m_pName, name);
  
   // SEt the value into the parameter
   m_value = initialValue;
      
   IncCtorCount();
} 

/******************************************************************************
IntParam::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void IntParamWorker::Write(FILE * pFile, int type) {

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%-13d  ", m_value);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%-13d  ", m_value);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s  ", m_pName);
      fprintf(pFile, "Value %d\n", m_value);

   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %d\n", m_pName, m_value);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-12s  ", m_pName);
   }/* end else() */
} /* end IntParam::WriteToFile() */
