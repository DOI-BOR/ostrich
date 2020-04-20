/******************************************************************************
File      : RespVarABC.h
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

Encapsulates an abstract response variable, the optimization analog of 
observations. This base class provides a unified interface between constraint 
and cost classes and regular and tied response variables.

Version History
01-10-05    lsm   created
******************************************************************************/
#ifndef RESP_VAR_ABC_H
#define RESP_VAR_ABC_H

#include "MyHeaderInc.h"

/******************************************************************************
class RespVarABC

 This class represents a generic response variable. A response variable is used
 in the computation of constraints and/or system cost.
******************************************************************************/
class RespVarABC
{
   public:
     virtual ~RespVarABC(void){ DBG_PRINT("RespVarABC::DTOR"); }
     virtual void Destroy(void) = 0;
     virtual void Write(FILE * pFile, int type) = 0;
     virtual double GetInitialVal(void) = 0;
     virtual double GetCurrentVal(void) = 0;
     virtual void SetCurrentVal(double val) = 0;
     virtual UnchangeableString GetName(void) = 0;     
}; /* end class RespVarABC */

#endif /* RespVarABC */
