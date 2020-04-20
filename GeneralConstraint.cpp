/******************************************************************************
File      : GeneralConstraint.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a single general constraint. General constraints are specified in 
the response variables group (ResponseVarGroup). The penalty is computed as the 
absolute value of the violation of the constraint multipled by a conversion factor 
which converts the units of the violation to a cost unit (dollars). That is, 
the conversion factor specifies the cost per unit of violation.

Version History
05-12-04    lsm   created
01-10-05    lsm   Class now interfaces with abstract resp. variables (RespVarABC)
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ConstraintABC.h"
#include "GenConstrainedOpt.h"
#include "RespVarABC.h"

#include "Exception.h"

/******************************************************************************
GetResponseVar()
******************************************************************************/
double GeneralConstraint::GetResponseVar(void)
{
	return m_pLoc->GetCurrentVal(); 
} /* end GetResponseVar() */

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
GeneralConstraint::GeneralConstraint
(IroncladString name, RespVarABC * pVar, double lwr, double upr, double conv)
{
   int len;
   
   len = (int)strlen(name) + 1;
   
   NEW_PRINT("char", len);
   m_Name = new char[len];

   NEW_PRINT("char", 40);
   m_TypeStr = new char[40];
   
   MEM_CHECK(m_TypeStr);

   strcpy(m_Name, name);
   strcpy(m_TypeStr, "General");

   m_pLoc = pVar;
   m_Upr = upr;
   m_Lwr= lwr;
   m_Conv = conv;
   m_Viol = 0.00;
   m_pNext = NULL;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Free up memory of the constraint and all other constraints in the list.
******************************************************************************/
void GeneralConstraint::Destroy(void)
{
   delete [] m_Name;
   delete [] m_TypeStr;
   delete m_pNext;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
AddConstraint()

Insert a constraint at the end of the linked list.
******************************************************************************/
void GeneralConstraint::AddConstraint(ConstraintABC * pNxt)
{ 
   if(m_pNext == NULL){ m_pNext = pNxt;} 
   else{ m_pNext->AddConstraint(pNxt);}
}/* end AddConstraint() */

/******************************************************************************
CalcPenalty()

Calculate the constraint violation and associated penalty.
******************************************************************************/
double GeneralConstraint::CalcPenalty(void)
{
   double tmp;
   tmp = m_pLoc->GetCurrentVal();
   
   m_Viol = 0.00;
   if(tmp < m_Lwr){ m_Viol = fabs(tmp-m_Lwr); }
   if(tmp > m_Upr){ m_Viol = fabs(tmp-m_Upr); }

   return (m_Viol*m_Conv);
} /* end CalcPenalty() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void GeneralConstraint::Write(FILE * pFile, int type)
{
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         fprintf(pFile, "%-12s  %E  %E  ", m_Name, m_Viol, m_Viol*m_Conv);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
	      fprintf(pFile, "%-12s  %.6lf  %.6lf  ", m_Name, m_Viol, m_Viol*m_Conv);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_BNR) :
      {
	      fprintf(pFile, "Name           Violation      Penalty        ");
         break;
      }/* end case(WRITE_BNR) */
      default:
      case(WRITE_DBG) :
      {
         fprintf(pFile, "******Constraint******\n");
         fprintf(pFile, "Name       : %s\n", m_Name);
         fprintf(pFile, "Type       : %s\n", m_TypeStr);
         fprintf(pFile, "Lower      : %.6lf     Upper     : %.6lf\n", m_Lwr, m_Upr);
         fprintf(pFile, "Conversion : %.6lf     Violation : %.6lf\n", m_Conv, m_Viol);
         fprintf(pFile, "Penalty    : %.6lf\n", m_Viol*m_Conv);
         m_pLoc->Write(pFile, type);
         break;
      }/* end case(WRITE_SCI) */
   }/* end switch() */
} /* end Write() */

