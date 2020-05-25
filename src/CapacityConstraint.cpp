/******************************************************************************
File      : CapacityConstraint.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Capacity constraints limit the summed value of a group of input parameters.
(for example, limits may be placed on the total pumping rate to ensure that an 
existing treatment plant is not overloaded). Constraint variables are 
stored in the ParameterGroup list and are identified by the name list. The 
penalty is computed as the absolute value of the violation of the constraint 
multiplied by a conversion factor which converts the units of the capacity violation 
(e.g. Length^3/Time for pumping rate) to a cost unit (dollars). That is, the 
conversion factor specifies the  cost per unit of capacity violation.

Version History
05-12-04    lsm   created
01-10-05    lsm   Generalized the PatoConstraintABC class
******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ConstraintABC.h"
#include "GenConstrainedOpt.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ResponseVarGroup.h"

#include "Exception.h"

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
CapacityConstraint::CapacityConstraint
(IroncladString name, IroncladString * nameList, int numNames, 
 ParameterGroup * pGroup, double lwr, double upr, double conv)
{
   char msg[DEF_STR_SZ];
   int i, len;

   len = (int)strlen(name) + 1;
   NEW_PRINT("char", len);
   m_Name = new char[len];

   NEW_PRINT("char", 50);
   m_TypeStr = new char[50];      
   MEM_CHECK(m_TypeStr);

   strcpy(m_Name, name);
   strcpy(m_TypeStr, "Capacity");
   
   m_NumVars = numNames;
   NEW_PRINT("ParameterABC *", m_NumVars);
   m_pParams = new ParameterABC *[m_NumVars];   
   MEM_CHECK(m_pParams);

   //determine indices by examining parameter names
   for(i = 0; i < m_NumVars; i++)
   {
      m_pParams[i] = pGroup->GetParamPtr(nameList[i]);

      //if name not found in list, must abort program
      if(m_pParams[i] == NULL)
      {
         sprintf(msg, "CapacityConstraint, unknown parameter : |%s|", nameList[i]);
         LogError(ERR_FILE_IO, msg);
         ExitProgram(1);         
      }/* end if() */
   }/* end for() */

   m_Upr = upr;
   m_Lwr= lwr;
   m_Conv = conv;
   m_Viol = 0.00;
   m_pNext = NULL;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
AddConstraint()

Insert a constraint at the end of the linked list.
******************************************************************************/
void CapacityConstraint::AddConstraint(ConstraintABC * pNxt)
{ 
   if(m_pNext == NULL){ m_pNext = pNxt;} 
   else{ m_pNext->AddConstraint(pNxt);}
}/* end AddConstraint() */

/******************************************************************************
Destroy()

Free up memory of the constraint and all other constraints in the list.
******************************************************************************/
void CapacityConstraint::Destroy(void)
{
   delete [] m_Name;
   delete [] m_TypeStr;
   delete [] m_pParams;
   delete m_pNext;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CalcPenalty()

Calculate the constraint violation and associated penalty.
******************************************************************************/
double CapacityConstraint::CalcPenalty(void)
{
   double cur, total;
   int i;

   total = 0.00;
   for(i = 0; i < m_NumVars; i++)
   {
      cur = m_pParams[i]->GetEstVal();
      total += cur;
   }/* end for() */

   m_Viol = 0.00;
   if(total < m_Lwr){ m_Viol = fabs(total-m_Lwr); }
   if(total > m_Upr){ m_Viol = fabs(total-m_Upr); }

   return (fabs(m_Viol)*m_Conv);
} /* end CalcPenalty() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void CapacityConstraint::Write(FILE * pFile, int type)
{
   int i;
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

         for(i = 0; i < m_NumVars; i++){ m_pParams[i]->Write(pFile, type);}
         break;
      }/* end case(WRITE_SCI) */
   }/* end switch() */
} /* end Write() */

