/******************************************************************************
File      : HydGradConstraint.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a single hydraulic gradient constraint. Hydraulic gradient 
constraints are composed of two head values which are stored in the response 
variables group (ResponseVarGroup). The difference between these two heads is 
the hydraulic gradient, which must be greater than or less than some constraint 
value. The penalty is computed as the absolute value of the violation 
of the constraint multipled by a conversion factor which converts the units of 
the gradient violation (Length) to a cost unit (dollars). That is, the conversion
factor specifies the cost per unit length of gradient violation.

Version History
05-12-04    lsm   created
01-10-05    lsm   Generalized the PatoConstraintABC class and modified to interface 
                  with abstract response variables (RespVarABC)
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
double HydGradConstraint::GetResponseVar(void)
{
	return m_pHead1->GetCurrentVal() - m_pHead2->GetCurrentVal(); 
}/* end GetResponseVar() */

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
HydGradConstraint::HydGradConstraint
(IroncladString name, RespVarABC * pHead1, RespVarABC * pHead2, double lwr, 
 double upr, double conv)
{
   int len;

   len = (int)strlen(name) + 1;

   NEW_PRINT("char", len);
   m_Name = new char[len];

   NEW_PRINT("char", 40);
   m_TypeStr = new char[40];

   MEM_CHECK(m_TypeStr);

   strcpy(m_Name, name);
   strcpy(m_TypeStr, "Hydraulic Gradient");

   m_pHead1 = pHead1;
   m_pHead2 = pHead2;
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
void HydGradConstraint::Destroy(void)
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
void HydGradConstraint::AddConstraint(ConstraintABC * pNxt)
{ 
   if(m_pNext == NULL){ m_pNext = pNxt;} 
   else{ m_pNext->AddConstraint(pNxt);}
}/* end AddConstraint() */

/******************************************************************************
CalcPenalty()

Calculate the constraint violation and associated penalty.
******************************************************************************/
double HydGradConstraint::CalcPenalty(void)
{
   double h1, h2, diff;
   h1 = m_pHead1->GetCurrentVal();
   h2 = m_pHead2->GetCurrentVal();
   diff = h1 - h2;

   m_Viol = 0.00;
   if(diff < m_Lwr){ m_Viol = fabs(diff-m_Lwr); }
   if(diff > m_Upr){ m_Viol = fabs(diff-m_Upr); }

   return (fabs(m_Viol)*m_Conv);
} /* end CalcPenalty() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void HydGradConstraint::Write(FILE * pFile, int type)
{
   switch(type)
   {
      case(WRITE_SCI) : 
      {
         fprintf(pFile, "%-12s  %E  %E  ", m_Name, m_Viol, fabs(m_Viol)*m_Conv);
         break;
      }/* end case(WRITE_SCI) */
      case(WRITE_DEC) :
      {
	      fprintf(pFile, "%-12s  %.6lf  %.6lf  ", m_Name, m_Viol, fabs(m_Viol)*m_Conv);
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
         fprintf(pFile, "Penalty    : %.6lf\n", fabs(m_Viol)*m_Conv);
         m_pHead1->Write(pFile, type);
         m_pHead2->Write(pFile, type);
         break;
      }/* end case(WRITE_SCI) */
   }/* end switch() */
} /* end Write() */

