/******************************************************************************
File      : TiedParamABC.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a 'tied' parameter. Tied parameters are variables in the model 
which are computed from the values of one or more model parameters. The ABC
for the tied parameters encapsulates the interface used by other Ostrich 
modules, allowing various specific tied parameter relationships (linear, 
exponential, etc.) to be implemented as needed with minimal code change (just 
need to add the specific tied parameter class and some additional input file
parsing).

These specific tied-parameter classes are supported:
TiedParamLin1   : linear function of one parameter
TiedParamLin2   : linear function of two parameters
TiedParamExp    : base exponential function of one parameter
TiedParamLog    : logarithmic function of one parameter
TiedDistXY      : distance between two (x,y) parameters
TiedParamSimpleRatio  : simple ratio of two parameters (ax + b) / (cy + d)
TiedParamComplexRatio  : complex ratio of three parameters 
(Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H) / (Ixyz + Jxy + Kxz + Lyz + Mx + Ny + Oz + P)
TiedParamConstant     : parameter is assigned a constant value
TiedParamWsum  : weighted sum of one or more parameters

Version History
07-07-04    lsm   Created
09-13-04    lsm   added distance (TiedDistXY)
01-01-07    lsm   added ratios (TiedParamSimpleRatio and TiedParamComplexRatio)
******************************************************************************/
#ifndef TIED_PARAM_ABC_H
#define TIED_PARAM_ABC_H

#include "MyHeaderInc.h"

/******************************************************************************
class TiedParamABC

Abstract base class of a tied-parameter.
******************************************************************************/
class TiedParamABC
{
   public:      
      virtual ~TiedParamABC(void){ DBG_PRINT("TiedParamABC::DTOR"); }      
      virtual void Destroy(void) = 0;
      virtual void   GetValAsStr(UnmoveableString valStr) = 0;
      virtual void   Write(FILE * pFile, int type) = 0;
      virtual double GetEstVal(void) = 0;      
      virtual UnchangeableString GetName(void) = 0;      
}; /* end class TiedParamABC */

/******************************************************************************
class TiedParamLin1

Represents a linear function of one parameter (F = aX+b)
******************************************************************************/
class TiedParamLin1 : public TiedParamABC
{
   public:      
      TiedParamLin1(void);
      TiedParamLin1(IroncladString name, MetaParameter * p1, 
                    UnmoveableString configStr);
      ~TiedParamLin1(void){ DBG_PRINT("TiedParamLin1::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_Tie;
      double m_C1, m_C0; //coefficients      
}; /* end class TiedParamLin */

/******************************************************************************
class TiedParamLin2

Represents a linear function of two parameters (F = aX + bY + cXY + d)
******************************************************************************/
class TiedParamLin2 : public TiedParamABC
{
   public:      
      TiedParamLin2(void);
      TiedParamLin2(IroncladString name, MetaParameter * p1, MetaParameter * p2,
                    UnmoveableString configStr);
      ~TiedParamLin2(void){ DBG_PRINT("TiedParamLin2::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_Tie1;
      MetaParameter m_Tie2;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedParamLin2 */

/******************************************************************************
class TiedParamExp

Represents an exponential function of one parameter (F = a BASE^^(bX) + c)
******************************************************************************/
class TiedParamExp : public TiedParamABC
{
   public:      
      TiedParamExp(void);
      TiedParamExp(IroncladString name, MetaParameter * p1, 
                    UnmoveableString configStr);
      ~TiedParamExp(void){ DBG_PRINT("TiedParamExp::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_Tie;
      double m_Base;
      double m_C2, m_C1, m_C0; //coefficients      
}; /* end class TiedParamExp */

/******************************************************************************
class TiedParamLog

Represents a logarithmic function of one parameter (F = a LOG(bX+c) + d)
******************************************************************************/
class TiedParamLog : public TiedParamABC
{
   public:      
      TiedParamLog(void);
      TiedParamLog(IroncladString name, MetaParameter * p1, 
                    UnmoveableString configStr);
      ~TiedParamLog(void){ DBG_PRINT("TiedParamLog::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_Tie;
      double m_Base;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedParamLog */

/******************************************************************************
class TiedDistXY

Represents the distance between two (x,y) parameters.
******************************************************************************/
class TiedDistXY : public TiedParamABC
{
   public:      
      TiedDistXY(void);
      TiedDistXY(IroncladString name, MetaParameter * px1, MetaParameter * py1,
                                      MetaParameter * px2, MetaParameter * py2,
                                      UnmoveableString configStr);
      ~TiedDistXY(void){ DBG_PRINT("TiedDistXY::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_X1;
      MetaParameter m_Y1;
      MetaParameter m_X2;
      MetaParameter m_Y2;
}; /* end class TiedDistXY */

/******************************************************************************
class TiedParamSimpleRatio

Represents a ratio of linear functions of two parameters :
   F = (aX + b) / (cY + d)
******************************************************************************/
class TiedParamSimpleRatio : public TiedParamABC
{
   public:      
      TiedParamSimpleRatio(void);
      TiedParamSimpleRatio(IroncladString name, MetaParameter * p1, MetaParameter * p2,
                     UnmoveableString configStr);
      ~TiedParamSimpleRatio(void){ DBG_PRINT("TiedParamSimpleRatio::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_Tie1;
      MetaParameter m_Tie2;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedParamSimpleRatio */


/******************************************************************************
class TiedParamComplexRatio

Represents a complex ratio of linear functions of three parameters :
   (Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H) 
   -------------------------------------------
   (Ixyz + Jxy + Kxz + Lyz + Mx + Ny + Oz + P)
******************************************************************************/
class TiedParamComplexRatio : public TiedParamABC
{
   public:      
      TiedParamComplexRatio(void);
      TiedParamComplexRatio(IroncladString name, MetaParameter * p1, MetaParameter * p2,
                     MetaParameter * p3, UnmoveableString configStr);
      ~TiedParamComplexRatio(void){ DBG_PRINT("TiedParamComplexRatio::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter m_X;
      MetaParameter m_Y;
      MetaParameter m_Z;
      double m_N[8]; //numerator coefficients
      double m_D[8]; //denominator coefficients
}; /* end class TiedParamComplexRatio */

/******************************************************************************
class TiedParamConstant

Represents a constant parameter.
******************************************************************************/
class TiedParamConstant : public TiedParamABC
{
   public:      
      TiedParamConstant(void);
      TiedParamConstant(IroncladString name, UnmoveableString pVal);
      ~TiedParamConstant(void){ DBG_PRINT("TiedParamConstant::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      double m_val;
}; /* end class TiedParamConstant */

/******************************************************************************
class TiedParamWsum

Represents a weighted combination of parameters.
******************************************************************************/
class TiedParamWsum : public TiedParamABC
{
   public:      
      TiedParamWsum(void);
      TiedParamWsum(IroncladString name, MetaParameter * p1, int nun,
                    UnmoveableString configStr);
      ~TiedParamWsum(void){ DBG_PRINT("TiedParamWsum::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      MetaParameter * m_pTie;
      double * m_pWgt; //weights
      int m_Num; //number of entries
}; /* end class TiedParamWsum */

#endif /* TIED_PARAM_ABC_H */
