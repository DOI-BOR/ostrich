/******************************************************************************
File      : ConstraintABC.h
Author    : L. Shawn Matott 
Copyright : 2016, L. Shawn Matott

Version History
11-29-16    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef CONSTRAINT_ABC_H
#define CONSTRAINT_ABC_H

#include "MyHeaderInc.h"

//forward decs
class ParameterGroup;
class RespVarABC;
class ParameterABC;

/******************************************************************************
class ConstraintABC (optimization constraints) base class
******************************************************************************/
class ConstraintABC
{
   public :
      virtual ~ConstraintABC(void){ DBG_PRINT("ConstraintABC::DTOR"); }
      virtual void Destroy(void)=0;
      virtual double CalcPenalty(void)=0;      
      virtual ConstraintABC * GetNext(void)=0;
      virtual void AddConstraint(ConstraintABC * pNxt)=0;
	   virtual void Write(FILE * pFile, int type)=0; 
	   virtual double GetLowerLimit(void)=0; 
	   virtual double GetUpperLimit(void)=0;
	   virtual double GetResponseVar(void)=0;
	   virtual UnchangeableString GetName(void)=0;
}; /* end class ConstraintABC */

/******************************************************************************
class GeneralConstraint

General constraints are imposed directly on the value a response variable
specified in the response variables group (ResponseVarGroup). The penalty is
computed as the absolute value of the violation of the constraint multipled by a
conversion factor which converts the units of the constraint to a cost unit
(dollars). That is, the conversion factor specifies the cost per unit of
violation.
******************************************************************************/
class GeneralConstraint : public ConstraintABC
{
public:
	~GeneralConstraint(void){ DBG_PRINT("GeneralConstraint::DTOR"); Destroy(); }
	void Destroy(void);
	GeneralConstraint(IroncladString name, RespVarABC * pVar, double lwr,
		double upr, double conv);
	double CalcPenalty(void);
	ConstraintABC * GetNext(void){ return m_pNext; }
	void AddConstraint(ConstraintABC * pNxt);
	void Write(FILE * pFile, int type);
	double GetLowerLimit(void){ return m_Lwr; }
	double GetUpperLimit(void){ return m_Upr; }
	double GetResponseVar(void);
	UnchangeableString GetName(void){ return m_Name; }

private:
	ConstraintABC * m_pNext;
	StringType m_Name;
	StringType m_TypeStr;
	//pointer into the response variable group
	RespVarABC * m_pLoc;
	double m_Lwr; //lower bound of the constraint
	double m_Upr; //upper bound of the constraint
	double m_Conv; //conversion factor (cost per unit violation)
	double m_Viol; //constraint violation
}; /* end class GeneralConstraint */

/******************************************************************************
class CapacityConstraint (capacity constraint)

Capacity constraints limit the summed value of a group of input parameters.
(for example, limits may be placed on the total pumping rate to ensure that an
existing treatment plant is not overloaded). Constraint variables are
stored in the ParameterGroup list and are identified by the name list. The
penalty is computed as the absolute value of the violation of the constraint
multiplied by a conversion factor which converts the units of the capacity violation
(e.g. Length^3/Time for pumping rate) to a cost unit (dollars). That is, the
conversion factor specifies the  cost per unit of capacity violation.
******************************************************************************/
class CapacityConstraint : public ConstraintABC
{
public:
	~CapacityConstraint(void){ DBG_PRINT("CapacityConstraint::DTOR"); Destroy(); }
	void Destroy(void);
	CapacityConstraint(IroncladString name, IroncladString * nameList,
		int numNames, ParameterGroup * pGroup, double lwr,
		double upr, double conv);
	double CalcPenalty(void);
	ConstraintABC * GetNext(void){ return m_pNext; }
	void AddConstraint(ConstraintABC * pNxt);
	void Write(FILE * pFile, int type);
	double GetLowerLimit(void){ return m_Lwr; }
	double GetUpperLimit(void){ return m_Upr; }
	double GetResponseVar(void){ return 0.00; }
	UnchangeableString GetName(void){ return m_Name; }

private:
	ConstraintABC * m_pNext;
	StringType m_Name;
	StringType m_TypeStr;
	ParameterABC ** m_pParams; //array of parameters in the ParameterGroup
	int m_NumVars;  //number of parameters in capacity summation
	double m_Lwr; //lower bound of the constraint
	double m_Upr; //upper bound of the constraint
	double m_Conv; //conversion factor (cost per unit violation)
	double m_Viol; //constraint violation
}; /* end class CapacityConstraint */

/******************************************************************************
class HydGradConstraint (hydraulic gradient constraint)

Hydraulic gradient constraints are composed of two head values which are stored
in the response variables group (ResponseVarGroup). The difference between these
two heads is the hydraulic gradient, which must be greater than or less than some
constraint value. The penalty is computed as the absolute value of the violation
of the constraint multipled by a conversion factor which converts the units of
the gradient violation (Length) to a cost unit (dollars). That is, the conversion
factor specifies the cost per unit length of gradient violation.
******************************************************************************/
class HydGradConstraint : public ConstraintABC
{
public:
	~HydGradConstraint(void){ DBG_PRINT("HydGradConstraint::DTOR"); Destroy(); }
	HydGradConstraint(IroncladString name, RespVarABC * pHead1,
		RespVarABC * pHead2, double lwr, double upr, double conv);
	void Destroy(void);
	double CalcPenalty(void);
	ConstraintABC * GetNext(void){ return m_pNext; }
	void AddConstraint(ConstraintABC * pNxt);
	void Write(FILE * pFile, int type);
	double GetLowerLimit(void){ return m_Lwr; }
	double GetUpperLimit(void){ return m_Upr; }
	double GetResponseVar(void);
	UnchangeableString GetName(void){ return m_Name; }

private:
	ConstraintABC * m_pNext;
	StringType m_Name;
	StringType m_TypeStr;
	//pointers into the response variable group
	RespVarABC * m_pHead1;
	RespVarABC * m_pHead2;
	double m_Lwr; //lower bound of the constraint
	double m_Upr; //upper bound of the constraint
	double m_Conv; //conversion factor (cost per unit violation)
	double m_Viol; //constraint violation
}; /* end class HydGradConstraint */

/******************************************************************************
class DrawdownConstraint (drawdown constraint)

Drawdown constraints are composed of the initial and current head values and are
enforced at user-specified locations as specified in the response variables group
(ResponseVarGroup). The difference between the initial and current heads is the
drawdown, which must be greater than or less than some constraint value. The
penalty is computed as the absolute value of the violation of the constraint
multipled by a conversion factor which converts the units of the drawdown violation
(Length) to a cost unit (dollars). That is, the conversion factor specifies the
cost per unit length of drawdown violation.
******************************************************************************/
class DrawdownConstraint : public ConstraintABC
{
public:
	~DrawdownConstraint(void){ DBG_PRINT("DrawdownConstraint::DTOR"); Destroy(); }
	void Destroy(void);
	DrawdownConstraint(IroncladString name, RespVarABC * pLoc, double lwr,
		double upr, double conv);
	double CalcPenalty(void);
	ConstraintABC * GetNext(void){ return m_pNext; }
	void AddConstraint(ConstraintABC * pNxt);
	void Write(FILE * pFile, int type);
	double GetLowerLimit(void){ return m_Lwr; }
	double GetUpperLimit(void){ return m_Upr; }
	double GetResponseVar(void);
	UnchangeableString GetName(void){ return m_Name; }

private:
	ConstraintABC * m_pNext;
	StringType m_Name;
	StringType m_TypeStr;
	//pointer into the response variable group where the drawdown location is stored
	RespVarABC * m_pLoc;
	double m_Lwr; //lower bound of the constraint
	double m_Upr; //upper bound of the constraint
	double m_Conv; //conversion factor (cost per unit violation)
	double m_Viol; //constraint violation
}; /* end class DrawdownConstraint */

/******************************************************************************
class ParticleCaptureConstraint (particle capture constraint)

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
******************************************************************************/
class ParticleCaptureConstraint : public ConstraintABC
{
public:
	~ParticleCaptureConstraint(void){ DBG_PRINT("ParticleCaptureConstraint::DTOR"); Destroy(); }
	void Destroy(void);
	ParticleCaptureConstraint(IroncladString name, RespVarABC * pX,
		RespVarABC * pY, Point2D * pPlume,
		int nv, double conv);
	double CalcPenalty(void);
	ConstraintABC * GetNext(void){ return m_pNext; }
	void AddConstraint(ConstraintABC * pNxt);
	void Write(FILE * pFile, int type);
	double GetLowerLimit(void){ return NEARLY_ZERO; }
	double GetUpperLimit(void){ return NEARLY_HUGE; }
	double GetResponseVar(void){ return 0.00; }
	UnchangeableString GetName(void){ return m_Name; }

private:
	ConstraintABC * m_pNext;
	StringType m_Name;
	StringType m_TypeStr;
	//pointers to the x- and y-coordinates in the response variable group
	RespVarABC * m_pXcoord;
	RespVarABC * m_pYcoord;
	Point2D * m_pPlume; //vertices in the plume geometry
	int m_NumVert;  //number of plume vertices
	double m_Conv; //conversion factor (cost per unit violation)
	double m_Viol; //constraint violation
}; /* end class ParticleCaptureConstraint */

#endif /* CONSTRAINT_ABC_H */






