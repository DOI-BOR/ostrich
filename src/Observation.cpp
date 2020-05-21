/******************************************************************************
File      : Observation.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates a single observation point.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
07-08-04    lsm   added break statements to all switch() cases
08-17-04    lsm   added reporting of memory allocations
01-01-07    lsm   added copy CTOR and Reconfigure() routines to support Surrogate-
                  model approach
******************************************************************************/
#include <stdio.h>
#include <string.h>

#include "Observation.h"
#include "ObjectiveFunction.h"

#include "Exception.h"

/******************************************************************************
CalcResidual()
  
Returns the possibly transformed and weighted residual at the observation point.
******************************************************************************/
double Observation::CalcResidual(bool bTransformed, bool bWeighted)
{
  double x, y;
  if(bTransformed == true) //transformed implies also weighted
  {
    x = BoxCox(m_MeasuredVal*m_Weight);
    y = BoxCox(m_ComputedVal*m_Weight);
  }
  else if(bWeighted == true)
  {
    x = (m_MeasuredVal*m_Weight);
    y = (m_ComputedVal*m_Weight);
  }
  else
  {
    x = (m_MeasuredVal);
    y = (m_ComputedVal);
  }
  return (x - y);
} /* end CalcResidual() */

/******************************************************************************
GetName()

Returns the name of the observation.
******************************************************************************/
UnchangeableString Observation::GetName(void)
{
  return m_Name;
} /* end GetName() */

/******************************************************************************
GetGroup()

Returns the group name of the observation.
******************************************************************************/
UnchangeableString Observation::GetGroup(void)
{
  return m_Group;
} /* end GetGroup() */

/******************************************************************************
CTOR

Dummy constructor of the class.
******************************************************************************/
Observation::Observation(void)
{
   m_Name = NULL;
   m_FileName = NULL;
   m_Keyword = NULL;
   m_Group = NULL;
   m_Tok = ' ';
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()
******************************************************************************/
void Observation::Destroy(void)
{
   if(m_Name != NULL)
   {
      delete [] m_Name;
      delete [] m_FileName;
      delete [] m_Keyword;
      delete [] m_Group;
   } /* end if() */

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
SetComputedVal()

Sets the computed value to the given value.
******************************************************************************/
void Observation::SetComputedVal(double computedVal)
{
  m_ComputedVal = computedVal;
} /* end SetComputedVal() */

/******************************************************************************
GetFileName()

Returns the file name associated with the observation. 
******************************************************************************/
UnchangeableString Observation::GetFileName(void)
{
  return m_FileName;
} /* end GetFileName() */

/******************************************************************************
GetKeyword()

Returns the key word associated with the observation. The extraction of the 
observation value depends on the key word as the extracting method first looks 
for the keyword.
******************************************************************************/
UnchangeableString Observation::GetKeyword(void)
{
  return m_Keyword;
} /* end GetKeyword() */

/******************************************************************************
GetLine()

Returns the line number associated with the observation.
******************************************************************************/
int Observation::GetLine(void)
{
  return m_Line;
} /* end GetLine() */

/******************************************************************************
GetColumn()

Returns the column  associated with the observation
******************************************************************************/
int Observation::GetColumn(void)
{
  return m_Column;
} /* end GetColumn() */

/******************************************************************************
GetMeasuredVal()

Returns the measured (observed) value for the observation point.
******************************************************************************/
double Observation::GetMeasuredVal(bool bTransformed, bool bWeighted)
{
  double x;
  if(bTransformed == true) //transformed implies also weighted
  {
    x = BoxCox(m_MeasuredVal*m_Weight);
  }
  else if(bWeighted == true)
  {
    x = (m_MeasuredVal*m_Weight);
  }
  else
  {
    x = (m_MeasuredVal);
  }
  return (x);
} /* end GetMeasuredVal() */

/******************************************************************************
GetComputedVal()

Returns the model computed value for the observation
******************************************************************************/
double Observation::GetComputedVal(bool bTransformed, bool bWeighted)
{
  double y;
  if(bTransformed == true) //transformed implies also weighted
  {
    y = BoxCox(m_ComputedVal*m_Weight);
  }
  else if(bWeighted == true)
  {
    y = (m_ComputedVal*m_Weight);
  }
  else
  {
    y = (m_ComputedVal);
  }
  return (y);
} /* end GetComputedValue() */


/******************************************************************************
GetObsWeight() 

Retrieve observation weight.
******************************************************************************/
double GetObsWeight(Observation * pObs)
{
   return pObs->GetWeight();
}/* end GetObsWeight() */

/******************************************************************************
CTOR 

Constructor for the class
******************************************************************************/
Observation::Observation
(
   IroncladString name,
   double measuredVal,
   double weight, 
	IroncladString fileName, 
   IroncladString keyword, 
   int line,
	int column,
   char tok,
   bool bAug,
   IroncladString group
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

   if(group != NULL)
   {
      len = (int)strlen(group) + 1;
      NEW_PRINT("char", len);
      m_Group = new char[len];
      strcpy(m_Group, group);
   }/* end if() */
   else
   {
      m_Group = new char[5];
      strcpy(m_Group, "none");
   }

   strcpy(m_Name, name);
   strcpy(m_Keyword, keyword);
   strcpy(m_FileName, fileName);

   m_MeasuredVal = measuredVal;   
   m_Weight      = weight;
   m_Line        = line;
   m_Column      = column;
   m_ComputedVal = 0.00;
   m_Tok = tok;
   m_bAug = bAug;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Copy CTOR 

Constructor for the class
******************************************************************************/
Observation::Observation(Observation * pCopy)
{   
   int len;

   len = (int)strlen(pCopy->GetName()) + 1;
   NEW_PRINT("char", len);
   m_Name = new char[len];   

   len = (int)strlen(pCopy->GetName()) + 1;
   NEW_PRINT("char", len);
   m_Keyword = new char[len];

   len = (int)strlen(pCopy->GetGroup()) + 1;
   NEW_PRINT("char", len);
   m_Group = new char[len];   

   len = (int)strlen(OST_OBS_FILE) + 1;
   NEW_PRINT("char", len);
   m_FileName = new char[len];

   MEM_CHECK(m_FileName);

   strcpy(m_Name, pCopy->GetName());
   strcpy(m_Keyword, pCopy->GetName());
   strcpy(m_FileName, OST_OBS_FILE);
   strcpy(m_Group, pCopy->GetGroup());

   m_MeasuredVal = pCopy->GetMeasuredVal(false, false);
   m_Weight      = pCopy->GetWeight();
   m_Line        = 0;
   m_Column      = 2;
   m_ComputedVal = 0.00;
   m_Tok = ' ';

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Reconfigure()

Reconfigure the parsing information.
******************************************************************************/
void Observation::Reconfigure
(
	IroncladString fileName, 
   IroncladString keyword, 
   int line,
	int column,
   char tok,
   bool bAug,
   IroncladString group
)
{   
   int len;

   delete [] m_FileName;
   delete [] m_Keyword;
   delete [] m_Group;

   len = (int)strlen(fileName) + 1;
   NEW_PRINT("char", len);
   m_FileName = new char[len];   
   
   len = (int)strlen(keyword) + 1;
   NEW_PRINT("char", len);
   m_Keyword = new char[len];

   MEM_CHECK(m_Keyword);

   if(group != NULL)
   {
      len = (int)strlen(group) + 1;
      NEW_PRINT("char", len);
      m_Group = new char[len];
      strcpy(m_Group, group);
   }/* end if() */
   else
   {
      m_Group = new char[5];
      strcpy(m_Group, "none");
   }

   strcpy(m_FileName, fileName);
   strcpy(m_Keyword, keyword);

   m_Line        = line;
   m_Column      = column;
   m_Tok = tok;
   m_bAug = bAug;
} /* end Reconfigure() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void Observation::Write(FILE * pFile, int type)
{
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         fprintf(pFile, "%E  %E  ", m_MeasuredVal, m_ComputedVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
	      fprintf(pFile, "%.6lf  %.6lf  ", m_MeasuredVal, m_ComputedVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_BNR) :
      {
	      fprintf(pFile, "%-23s  measured       computed       ", m_Name);
         break;
      }/* end case(WRITE_SCI) */
      default:
      case(WRITE_DBG) :
      {
	      fprintf(pFile, "%s  %E  %E  %s  %s  %d  %d  %c %E %s\n",
                 m_Name, m_MeasuredVal, m_Weight, m_FileName, 
		           m_Keyword, m_Line, m_Column, m_Tok, m_ComputedVal, m_Group);
         break;
      }/* end case(WRITE_DBG) */
   }/* end switch() */
} /* end Write() */

/******************************************************************************
WriteSim()

Writes simulated output to pFile.
******************************************************************************/
void Observation::WriteSim(FILE * pFile, int type)
{
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         fprintf(pFile, "%E  ", m_ComputedVal);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
	      fprintf(pFile, "%.6lf  ", m_ComputedVal);
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
	      fprintf(pFile, "%s  %E  %E  %s  %s  %d  %d  %c %E %s\n",
                 m_Name, m_MeasuredVal, m_Weight, m_FileName, 
		           m_Keyword, m_Line, m_Column, m_Tok, m_ComputedVal, m_Group);
         break;
      }/* end case(WRITE_DBG) */
   }/* end switch() */
} /* end WriteSim() */
