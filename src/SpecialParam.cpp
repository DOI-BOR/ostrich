/******************************************************************************
File      : SpecialParam.cpp
Author    : L. Shawn Matott
Copyright : 2009, L. Shawn Matott

Encapsulates special Ostrich parameters. Special parameters are variables 
the Ostrich computes internally which may be used to trigger model pre-emption.

For example, if the current best cost function is exceeded then it might make
sense to halt model early.

Version History
05-14-09    lsm   added copyright information and initial comments.
******************************************************************************/
#include <stdio.h>
#include <string.h>

#include "ParameterABC.h"
#include "ConstraintABC.h"
#include "GenConstrainedOpt.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
SpecialParam::GetValAsStr()
******************************************************************************/
void SpecialParam::GetValAsStr(UnmoveableString valStr)
{
	GetPreciseValAsStr(valStr, GetTransformedVal()); 
}/* end GetValAsStr() */

/******************************************************************************
SpecialParam::Destroy() 
******************************************************************************/
void SpecialParam::Destroy(void)
{
   delete [] m_pName;
   delete [] m_pType; 
   delete [] m_pLimit;
   delete [] m_pConstraint;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (SpecialParam)
******************************************************************************/
SpecialParam::SpecialParam(void)
{
   m_pName = NULL;
   m_pType = NULL; 
   m_pLimit = NULL;
   m_pConstraint = NULL;

   m_EstVal = -1.00;
   m_bSet = false;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
SpecialParam::SetEstVal()

Sets the estimated value of the parameter, depending on it's type (cost or constraint).
******************************************************************************/
void SpecialParam::SetEstVal(double minObj, double minCon)
{
   if(!m_bSet) return; //use initial values until algorithm is ready....

   if(strcmp(m_pType, "BestCost") == 0)
   {
	  m_EstVal = minObj;
   }
   else if(strcmp(m_pType, "BestConstraint") == 0)
   {
      double limit;
      if(strcmp(m_pLimit, "upper") == 0)
      {
         limit = GetConstraint()->GetUpperLimit();
         m_EstVal = MyMax(limit, minCon);
      }
      else if(strcmp(m_pLimit, "lower") == 0)
      {
         limit = GetConstraint()->GetLowerLimit();
         m_EstVal = MyMin(limit, minCon);
      }
      else
      {
         m_EstVal = minCon;
      }
   }
} /* end SpecialParam::SetEstVal() */

/******************************************************************************
SpecialParam::GetConstraint()
******************************************************************************/
ConstraintABC * SpecialParam::GetConstraint(void)
{
	if(strcmp(m_pType, "BestConstraint") == 0)
	{
		return(GetConstraintByName(m_pConstraint));
	}
	return NULL;
}/* end GetConstraint() */

/******************************************************************************
CTOR (SpecialParam)
******************************************************************************/
SpecialParam::SpecialParam(IroncladString name, IroncladString type, 
                           IroncladString limit, IroncladString constraint,
						   double init)
{
   int len;

   len = (int)strlen(name) + 10;  
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   len = (int)strlen(type) + 10;  
   NEW_PRINT("char", len);
   m_pType = new char[len];
   MEM_CHECK(m_pType);

   strcpy(m_pType, type);

   len = (int)strlen(limit) + 10;  
   NEW_PRINT("char", len);
   m_pLimit = new char[len];
   MEM_CHECK(m_pLimit);

   strcpy(m_pLimit, limit);

   len = (int)strlen(constraint) + 10;  
   NEW_PRINT("char", len);
   m_pConstraint = new char[len];
   MEM_CHECK(m_pConstraint);

   strcpy(m_pConstraint, constraint);

   m_EstVal = init;
   m_bSet = false;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
SpecialParam::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void SpecialParam::Write(FILE * pFile, int type)
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
      fprintf(pFile, "Name  = %s  ", m_pName);
	  fprintf(pFile, "Type  = %s  ", m_pType);
	  fprintf(pFile, "Limit = %s  ", m_pLimit);
	  fprintf(pFile, "Constraint = %s  ", m_pConstraint);
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
} /* end SpecialParam::Write() */
