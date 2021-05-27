/******************************************************************************
File      : ParameterABC.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a parameter. Parameters are variables in the model which are to 
be calibrated.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-25-03    lsm   Modified to support multiple stages of unit conversion.
07-05-04    lsm   added integer and combinatorial parameter support
12-02-04    lsm   made ConvertInVal() a required member function
03-03-05    lsm   Added support for ON/OFF parameter threshold
******************************************************************************/
#ifndef PARAMETER_WORKER
#define PARAMETER_WORKER

#include "MyHeaderInc.h"
#include "string"

// forward decs
class ConstraintABC;

/******************************************************************************
class ParameterABC

Abstract base class of a parameter.
******************************************************************************/
class ParameterWorkerABC
{
public:
    virtual ~ParameterWorkerABC(void) { DBG_PRINT("ParameterABC::DTOR"); }
    virtual void Destroy(void) = 0;

    virtual void Write(FILE* pFile, int type) = 0;
    virtual UnchangeableString GetName(void) = 0;
    virtual double GetValue(void) = 0;
    virtual void SetValue(double val) = 0;
    virtual void GetValAsStr(UnmoveableString valStr) = 0;
    virtual const char* GetType(void) = 0;
}; 


/******************************************************************************
class RealParamWorker

Represents a continuously varying parameter
******************************************************************************/
class RealParamWorker : public ParameterWorkerABC{
   public:      
       RealParamWorker(void);
       RealParamWorker(IroncladString name, double value, IroncladString fixFmt);     
      ~RealParamWorker(void){ DBG_PRINT("RealParamWorker::DTOR"); Destroy();}
      void Destroy(void);

      // Declare functions for values
      void   GetValAsStr(UnmoveableString valStr);
      UnchangeableString GetName(void) { return m_pName; }
      void   Write(FILE* pFile, int type);

      double GetValue(void){ return m_value;}
      void SetValue(double value) { m_value = value; }
      const char* GetType(void) { return "real"; };
      
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      double m_value;

}; 

/******************************************************************************
class IntParamWorker

Represents an integer parameter
******************************************************************************/
class IntParamWorker : public ParameterWorkerABC {
   public:      
       IntParamWorker(void);
       IntParamWorker(IroncladString name, int value);      
      ~IntParamWorker(void){ DBG_PRINT("IntParam::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){sprintf(valStr, "%d", m_value);}
      UnchangeableString GetName(void) { return m_pName; }
      void Write(FILE* pFile, int type);

      double GetValue(void){ return (double)m_value;}
      void SetValue(double value){ m_value = value; }
      const char * GetType(void) {return "integer";}


   private:
      StringType m_pName;
      int m_value;
      
}; 

#endif 
