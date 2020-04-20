/******************************************************************************
File      : DrawdownConstraint.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a single drawdown constraint. Drawdown constraints are composed of 
the initial and current head values and are enforced at user-specified locations 
as specified in the response variables group (ResponseVarGroup). The difference 
between the initial and current heads is the drawdown, which must be greater than 
or less than some constraint value. The penalty is computed as the absolute value 
of the violation of the constraint multipled by a conversion factor which converts 
the units of the drawdown violation  (Length) to a cost unit (dollars). That is, 
the conversion factor specifies the cost per unit length of drawdown violation.

Version History
05-12-04    lsm   created
01-10-05    lsm   Generalized the PatoConstraintABC class
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ConstraintABC.h"
#include "RespVarABC.h"

#include "Exception.h"

/******************************************************************************
GetResponseVar()
******************************************************************************/
double DrawdownConstraint::GetResponseVar(void)
{ 
	return m_pLoc->GetCurrentVal(); 
}/* end GetResponseVar() */

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
DrawdownConstraint::DrawdownConstraint
(IroncladString name, RespVarABC * pLoc, double lwr, double upr, double conv)
{
   int len;
   len = (int)strlen(name) + 1;
   NEW_PRINT("char", len);
   m_Name = new char[len];   
   NEW_PRINT("char", 40);
   m_TypeStr = new char[40];   
   MEM_CHECK(m_TypeStr);

   strcpy(m_Name, name);
   strcpy(m_TypeStr, "Drawdown");

   m_pLoc = pLoc;
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
void DrawdownConstraint::Destroy(void)
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
void DrawdownConstraint::AddConstraint(ConstraintABC * pNxt)
{ 
   if(m_pNext == NULL){ m_pNext = pNxt;} 
   else{ m_pNext->AddConstraint(pNxt);}
}/* end AddConstraint() */

/******************************************************************************
CalcPenalty()

Calculate the constraint violation and associated penalty.
******************************************************************************/
double DrawdownConstraint::CalcPenalty(void)
{
   double h1, h2, diff;
   h1 = m_pLoc->GetCurrentVal();
   h2 = m_pLoc->GetInitialVal();
   //drawdown is initial value minus current value (water level is decreaeing)
   diff = h2 - h1;

   m_Viol = 0.00;
   if(diff < m_Lwr){ m_Viol = fabs(diff-m_Lwr); }
   if(diff > m_Upr){ m_Viol = fabs(diff-m_Upr); }

   return (m_Viol*m_Conv);
} /* end CalcPenalty() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void DrawdownConstraint::Write(FILE * pFile, int type)
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

