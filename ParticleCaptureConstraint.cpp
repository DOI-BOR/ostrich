/******************************************************************************
File      : ParticleCaptureConstraint.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Particle capture constraints require that the location of a given particle be
within a well or within the original plume extents at the end of the planning 
horizon. These are specfied as (X,Y) pairs in the response variable group along
with a polygon that defines the plume geometry. At the end of the planning period,
a point-in-polygon test is performed to determine if the particle is in violation
of the capture/containment constraint. The penalty is computed as the square of 
the distance from the particle to the nearest plume boundary multiplied by a 
conversion factor which converts units from (Length^2) to cost (dollars).
Therefore, the conversion factor is the cost per unit violation of the particle
capture constraint.

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
#include "GeometryUtility.h"

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
ParticleCaptureConstraint::ParticleCaptureConstraint
(IroncladString name, RespVarABC * pX, RespVarABC * pY, Point2D * pPlume, 
 int nv, double conv)
{
   int len;

   len = (int)strlen(name) + 1;

   NEW_PRINT("char", len);
   m_Name = new char[len];

   NEW_PRINT("char", 50);
   m_TypeStr = new char[50];

   MEM_CHECK(m_TypeStr);

   strcpy(m_Name, name);
   strcpy(m_TypeStr, "Particle Capture");
   
   m_pXcoord = pX;
   m_pYcoord = pY;

   m_NumVert = nv;
   m_pPlume = pPlume;

   m_Conv = conv;
   m_Viol = 0.00;
   m_pNext = NULL;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Free up memory of the constraint and all other constraints in the list.
******************************************************************************/
void ParticleCaptureConstraint::Destroy(void)
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
void ParticleCaptureConstraint::AddConstraint(ConstraintABC * pNxt)
{ 
   if(m_pNext == NULL){ m_pNext = pNxt;} 
   else{ m_pNext->AddConstraint(pNxt);}
}/* end AddConstraint() */

/******************************************************************************
CalcPenalty()

Calculate the constraint violation and associated penalty.
******************************************************************************/
double ParticleCaptureConstraint::CalcPenalty(void)
{
   bool inPoly;
   Point2D loc;

   /*
   If the location of the particle is not inside the 
   plume, compute the square of the distance from the 
   particle to the nearest edge of the plume.
   */
   m_Viol = 0.00; //hoping that particle is in plume
   loc.x = m_pXcoord->GetCurrentVal();
   loc.y = m_pYcoord->GetCurrentVal();
   inPoly = PointInPoly(loc, m_pPlume, m_NumVert);
   if(inPoly == false){ m_Viol = DistToPoly(loc, m_pPlume, m_NumVert); }
   m_Viol *= m_Viol; //square result
   
   return (fabs(m_Viol)*m_Conv);
} /* end CalcPenalty() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void ParticleCaptureConstraint::Write(FILE * pFile, int type)
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
         fprintf(pFile, "Conversion : %.6lf     Violation : %.6lf\n", m_Conv, m_Viol);
         fprintf(pFile, "Penalty    : %.6lf\n", fabs(m_Viol)*m_Conv);
         fprintf(pFile, "------Plume Coords------\n");
         for(i = 0; i < m_NumVert; i++)
         {
            fprintf(pFile, "(%.6lf,%.6lf)\n", m_pPlume[i].x, m_pPlume[i].y);
         }/* end for() */

         m_pXcoord->Write(pFile, type);
         m_pYcoord->Write(pFile, type);
         break;
      }/* end case(WRITE_SCI) */
   }/* end switch() */
} /* end Write() */

