/******************************************************************************
File      : GenConstrainedOpt.h
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

Defines a general constrained optimization extension to the ObjectiveFunction class.

This class supports a veriety of cost and constraint formulations, allowing users to
define fairly generic objective functions without having to write a separate driver 
program.

This class instantiates a set of constraint classes which can be 
combined with the system cost using a user-selected penalty method (additive penalty, 
multiplicative penalty, etc.). Cost and constraints are made up of response variables
which are functions of model output and/or model parameters.
   
Version History
01-10-05    lsm   created
******************************************************************************/
#ifndef GCOP_H
#define GCOP_H

#include "MyHeaderInc.h"

// parent class
#include "ObjectiveFunction.h"

//forward declarations
class ConstraintABC;
class RespVarABC;
class ResponseVarGroup;
class ParameterGroup;

/******************************************************************************
class GCOP (General Constrained Optimization Problem)
******************************************************************************/
class GCOP : public ObjectiveFunction
{
   public :
      ~GCOP(void){ DBG_PRINT("GCOP::DTOR"); Destroy(); }
      void Destroy(void);
      GCOP(ParameterGroup * pParamGroup);      
      int CalcMultiObjFunc(double * pF, int nObj);
      double CalcObjFunc(void);
      void WriteSetupToFile(FILE * pFile);
      void WriteConstraints(FILE * pFile, int type);
      ConstraintABC * GetConstraintPtr(IroncladString pName);
	  void * GetResponseVarGroup(void);

   private :
      void InitFromFile(void);
      void InitResponseVars(void);
      void InitConstraints(void);

      LmtPenType  m_PenType;

      RespVarABC * m_pCostFunc;

      //multi-objective formulation
      RespVarABC ** m_pMultiObjCostFunc;
      int m_NumMultiObjCostFuncs;

      ParameterGroup * m_pParamGroup; // design variables
      ResponseVarGroup * m_pRespGroup; // response variables
      ConstraintABC * m_pConstraints; //linked-list of constraints
}; /* end class GCOP */

//C-style functions
extern "C"
{
   IroncladString GetPenMethStr(LmtPenType i);
}

#endif /* GCOP */
