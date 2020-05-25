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

#include "ParameterABC.h"

#include "Exception.h"
#include "FortranSupportUtilities.h"
#include "Utility.h"

/******************************************************************************
RealParam::Destroy() 
******************************************************************************/
void RealParam::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pFixFmt;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (RealParam)
******************************************************************************/
RealParam::RealParam(void)
{
   m_pName = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
RealParam::SetEstVal()

Sets the estimated value of the parameter. If the resquested estVal exceeds the
parameter bounds, the amount of violation is returned.
******************************************************************************/
double RealParam::SetEstVal(double estVal)
{
   double viol = 0.00;
   char msg[80];

   if(estVal < m_LwrBnd) 
   {            
      sprintf(msg, "%E < lower bound (%E)", estVal, m_LwrBnd);      
      LogError(ERR_PRM_BNDS, msg);
	  viol = m_LwrBnd - estVal;
      estVal = m_LwrBnd;
   }/* end if() */
   if(estVal > m_UprBnd) 
   {      
      sprintf(msg, "%E > upper bound (%E)", estVal, m_UprBnd);
      LogError(ERR_PRM_BNDS, msg);
	  viol = estVal - m_UprBnd;
      estVal = m_UprBnd;
   }/* end if() */

   //handle parameter threshold
   if((estVal < m_ThreshUpr) && (estVal > m_ThreshLwr))
   {
      estVal = m_ThreshOff;
   }

   m_EstVal = estVal;
   return viol;
} /* end RealParam::SetEstVal() */

/******************************************************************************
RealParam::GetTransformedVal()

Retrieves the current estimated value of the parameter, ready to be submitted
to the model executable, after transforming from log units, if nescessary.
******************************************************************************/
double RealParam::GetTransformedVal(void)
{
   switch(m_TransID[TX_OST])
   {
      case(TX_NONE)  : return m_EstVal;
      case(TX_LOG10) : return pow(10, m_EstVal);
      case(TX_LN)    : return exp(m_EstVal);
      default        : return 0.00;
   }/* end switch () */
} /* end RealParam::GetTransformedVal() */

/******************************************************************************
RealParam::ConvertOutVal()

Converts the input value based on users choice of output style.
******************************************************************************/
double RealParam::ConvertOutVal(double val)
{
   switch(m_TransID[TX_OST])
   {
      case(TX_NONE)  : 
         switch(m_TransID[TX_OUT])
         {
            case(TX_NONE)  : return val;
            case(TX_LOG10) : return log10(val);
            case(TX_LN)    : return log(val);
            default        : return 0.00;
         }/* end switch () */
      case(TX_LOG10) :          
         switch(m_TransID[TX_OUT])
         {
            case(TX_NONE)  : return pow(10, val);
            case(TX_LOG10) : return val;
            case(TX_LN)    : return log(pow(10,val));
            default        : return 0.00;
         }/* end switch () */
      case(TX_LN)    : 
         switch(m_TransID[TX_OUT])
         {
            case(TX_NONE)  : return exp(val);
            case(TX_LOG10) : return log10(exp(val));
            case(TX_LN)    : return val;
            default        : return 0.00;
         }/* end switch () */
      default        : return 0.00;
   }/* end switch () */
} /* end RealParam::ConvertOutVal() */

/******************************************************************************
RealParam::ConvertInVal()

Converts the input value based on users choice of Ostrich style. The resulting
value will be consistent with the transformation used by the Ostrich program.
******************************************************************************/
double RealParam::ConvertInVal(double val)
{
   switch(m_TransID[TX_OST])
   {
      case(TX_NONE)  : 
         switch(m_TransID[TX_IN])
         {
            case(TX_NONE)  : return val;
            case(TX_LOG10) : return pow(10, val);
            case(TX_LN)    : return exp(val);
            default        : return 0.00;
         }/* end switch () */
      case(TX_LOG10) :      
         switch(m_TransID[TX_IN])
         {
            case(TX_NONE)  : 
            {
                if (val <= 0.00) { val = NEARLY_ZERO;}
                return log10(val);
            }
            case(TX_LOG10) : return val;
            case(TX_LN)    : return log10(exp(val));
            default        : return 0.00;
         }/* end switch () */
      case(TX_LN)    : 
         switch(m_TransID[TX_IN])
         {
            case(TX_NONE)  : 
            {
               if (val <= 0.00) { val = NEARLY_ZERO;}
               return log(val);
            }              
            case(TX_LOG10) : return log(pow(10,val));
            case(TX_LN)    : return val;
            default        : return 0.00;
         }/* end switch () */
      default        : return 0.00;
   }/* end switch () */
} /* end RealParam::ConvertInVal() */

/******************************************************************************
CTOR (RealParam)
******************************************************************************/
RealParam::RealParam(IroncladString name, double initialValue, 
                     double lowerBound, double upperBound, IroncladString txIn, 
                     IroncladString txOst, IroncladString txOut, IroncladString fixFmt)
{
   int len;

   len = (int)strlen(name) + 10;  
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);
  
   SetTransformation(TX_IN, txIn);
   SetTransformation(TX_OST, txOst);
   SetTransformation(TX_OUT, txOut);

   m_InitVal = ConvertInVal(initialValue);
   m_ThreshLwr = m_ThreshUpr = m_ThreshOff = m_LwrBnd = ConvertInVal(lowerBound);
   m_UprBnd = ConvertInVal(upperBound);
   m_EstVal = ConvertInVal(initialValue);

   len = (int)strlen(fixFmt) + 10;  
   NEW_PRINT("char", len);
   m_pFixFmt = new char[len];
   MEM_CHECK(m_pFixFmt);

   strcpy(m_pFixFmt, fixFmt);
      
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
RealParam::SetTransformation()

Sets the ID of the transformation to be performed on the parameter prior to 
writing model input files.
******************************************************************************/
void RealParam::SetTransformation(TransformStageEnum which, 
                                  IroncladString transformation)
{   
   char errMsg[DEF_STR_SZ];

   m_TransID[which] = TX_NONE;
   if(strcmp(transformation, "none") == 0) {m_TransID[which] = TX_NONE;}
   else if(strcmp(transformation, "log10") == 0)  {m_TransID[which] = TX_LOG10;}
   else if(strcmp(transformation, "ln") == 0) {m_TransID[which] = TX_LN;}
   else //unknown
   {
      sprintf(errMsg, "Unknown transformation: %s", transformation);
      LogError(ERR_FILE_IO, errMsg);
   }
} /* end RealParam::SetTransformation() */

/******************************************************************************
RealParam::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void RealParam::Write(FILE * pFile, int type)
{
   double val;

   val = ConvertOutVal(m_EstVal);

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
      fprintf(pFile, "Name = %s  ", m_pName);
      fprintf(pFile, "Transformation= %d\n", m_TransID[TX_OST]);
      fprintf(pFile, "Initial Value %E\n", m_InitVal);
      fprintf(pFile, "Lower Bound %E\n", m_LwrBnd);
      fprintf(pFile, "Upper Bound %E\n", m_UprBnd);
      fprintf(pFile, "Lower Threshold %E\n", m_ThreshLwr);
      fprintf(pFile, "Upper Threshold %E\n", m_ThreshUpr);
      fprintf(pFile, "Off Threshold %E\n", m_ThreshOff);
      fprintf(pFile, "Est Value = %E\n", m_EstVal);
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
} /* end RealParam::WriteToFile() */

/******************************************************************************
RealParam::GetValAsStr()
******************************************************************************/
void RealParam::GetValAsStr(UnmoveableString valStr)
{
   bool bOk;
   double val = GetTransformedVal();
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
}/* end RealParam::GetValAsStr() */

/******************************************************************************
IntParam::Destroy()
******************************************************************************/
void IntParam::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (IntParam)
******************************************************************************/
IntParam::IntParam(void)
{
   m_pName = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
IntParam::SetEstVal()

Sets the estimated value of the parameter. If the resquested estVal exceeds the
parameter bounds, the amount of violation is returned.
******************************************************************************/
double IntParam::SetEstVal(double estVal)
{
   double viol = 0.00;
   int tmp;
   char msg[80];

   tmp = (int)(estVal + 0.5); //round up real value

   if(tmp < m_LwrBnd) 
   {            
      sprintf(msg, "%d < lower bound (%d)", tmp, m_LwrBnd);      
      LogError(ERR_PRM_BNDS, msg);
	  viol = (double)(m_LwrBnd - tmp);
      tmp = m_LwrBnd;
   }/* end if() */
   if(tmp > m_UprBnd) 
   {      
      sprintf(msg, "%d > upper bound (%d)", tmp, m_UprBnd);
      LogError(ERR_PRM_BNDS, msg);
	  viol = (double)(tmp - m_UprBnd);
      tmp = m_UprBnd;
   }/* end if() */

   //handle parameter threshold
   if((tmp < m_ThreshUpr) && (tmp > m_ThreshLwr))
   {
      tmp = m_ThreshOff;
   }

   m_EstVal = tmp;
   return viol;
} /* end IntParam::SetEstVal() */

/******************************************************************************
CTOR (IntParam)
******************************************************************************/
IntParam::IntParam(IroncladString name, int initialValue, int lowerBound, 
                   int upperBound)
{
   int len;

   len = (int)strlen(name) + 10;  
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);
  
   m_InitVal = initialValue;
   m_ThreshLwr = m_ThreshUpr = m_ThreshOff = m_LwrBnd = lowerBound;
   m_UprBnd = upperBound;
   m_EstVal = initialValue;
      
   IncCtorCount();
} /* end CTOR (IntParam) */

/******************************************************************************
IntParam::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void IntParam::Write(FILE * pFile, int type)
{
   int val;

   val = m_EstVal;

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%-13d  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%-13d  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s  ", m_pName);
      fprintf(pFile, "Initial Value %d\n", m_InitVal);
      fprintf(pFile, "Lower Bound  %d\n", m_LwrBnd);
      fprintf(pFile, "Upper Bound  %d\n", m_UprBnd);
      fprintf(pFile, "Lower Threshold %d\n", m_ThreshLwr);
      fprintf(pFile, "Upper Threshold %d\n", m_ThreshUpr);
      fprintf(pFile, "Threshold Off %d\n", m_ThreshOff);
      fprintf(pFile, "Est Value =  %d\n", m_EstVal);
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
} /* end IntParam::WriteToFile() */
