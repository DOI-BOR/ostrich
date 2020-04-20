/******************************************************************************
File      : TiedRespVar.h
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

Contains defintions for 'tied' response variables. Tied response variables 
are computed as functions of one or more response variables, which are read 
from model input and/or output files. The ABC for response variable 
encapsulates the interface used by other Ostrich modules, allowing various 
specific tied response variable relationships (linear, exponential, etc.) to 
be implemented as needed with minimal code change (just need to add the 
specific tied response variable class and some additional input file parsing).

These specific tied-response variable classes are supported:
TiedRespVarLin1  : linear function of one response variable
TiedRespVarLin2  : linear function of two response variables
TiedRespVarWsum  : weighted sum of one or more response variables

Version History
01-10-05    lsm   Created
01-20-05    lsm   added support for weighted sum tied resp. vars.
******************************************************************************/
#ifndef TIED_RESP_VAR_H
#define TIED_RESP_VAR_H

#include "MyHeaderInc.h"

// parent class
#include "RespVarABC.h"

/******************************************************************************
class TiedRespVarLin1

Represents a linear function of one resp. var. (F = aX+b)
******************************************************************************/
class TiedRespVarLin1 : public RespVarABC
{
   public:      
      TiedRespVarLin1(void);
      TiedRespVarLin1(IroncladString name, RespVarABC * p1, 
                    UnmoveableString configStr);
      ~TiedRespVarLin1(void){ DBG_PRINT("TiedRespVarLin1::DTOR"); Destroy();}
      void Destroy(void);

      void   Write(FILE * pFile, int type);      
      double GetCurrentVal(void);     
      double GetInitialVal(void);
      void SetCurrentVal(double curVal){ return;} //do nothing...
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      RespVarABC * m_pTie;
      double m_C1, m_C0; //coefficients      
}; /* end class TiedRespVarLin1 */

/******************************************************************************
class TiedRespVarLin2

Represents a linear function of two resp. vars (F = aX + bY + cXY + d)
******************************************************************************/
class TiedRespVarLin2 : public RespVarABC
{
   public:      
      TiedRespVarLin2(void);
      TiedRespVarLin2(IroncladString name, RespVarABC * p1, RespVarABC * p2,
                    UnmoveableString configStr);
      ~TiedRespVarLin2(void){ DBG_PRINT("TiedRespVarLin2::DTOR"); Destroy();}
      void Destroy(void);

      void   Write(FILE * pFile, int type);
      double GetCurrentVal(void);     
      double GetInitialVal(void);
      void SetCurrentVal(double curVal){ return;} //do nothing...
      UnchangeableString GetName(void){ return m_pName;}

   private:
      StringType m_pName;
      RespVarABC * m_pTie1;
      RespVarABC * m_pTie2;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedRespVarLin2 */

/******************************************************************************
class TiedRespVarWsum

Represents a weighted sum of resp. vars (F = w1*X1 + w2*X2 + ... + wn*Xn)
******************************************************************************/
class TiedRespVarWsum : public RespVarABC
{
   public:      
      TiedRespVarWsum(void);
      TiedRespVarWsum(IroncladString name, RespVarABC ** pList, int nrv, 
                      UnmoveableString configStr);
      ~TiedRespVarWsum(void){ DBG_PRINT("TiedRespVarWsum::DTOR"); Destroy();}
      void Destroy(void);

      void   Write(FILE * pFile, int type);
      double GetCurrentVal(void);     
      double GetInitialVal(void);
      void SetCurrentVal(double curVal){ return;} //do nothing...
      UnchangeableString GetName(void){ return m_pName;}

   private:
      StringType m_pName;
      RespVarABC ** m_pList; //response variables
      double * m_pWgt; //weights
      int m_Num;
}; /* end class TiedRespVarWsum */

#endif /* TIED_RESP_VAR_H */
