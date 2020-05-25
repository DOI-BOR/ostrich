/******************************************************************************
File      : ResponseVar.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a single response variable, the optimization analog of observations.

Version History
05-10-04    lsm   created
07-08-04    lsm   added support for dynamic row and column determination
01-10-05    lsm   Modified to inherit from abstract response variables (RespVarABC)
******************************************************************************/
#ifndef RESPONSE_VAR_H
#define RESPONSE_VAR_H

#include "MyHeaderInc.h"

// parent class
#include "RespVarABC.h"

// forward decs
class ParameterABC;
class TiedParamABC;

/******************************************************************************
class ResponseVar

 This class represents a response variable. A response variable is read from the
 output files of the model executable and is used in computation of the objective
 function.

 Each response variable has parsing paramters: fileName, keyword, line, and 
 column such that the  program associates each response variable with the value 
 which is found on the line and column after the first occurence of the keyword 
 in <fileName>.
******************************************************************************/
class ResponseVar : public RespVarABC
{
   public:
      ResponseVar(IroncladString name, IroncladString fileName, 
                  IroncladString keyword, int line, int column, char tok, bool bAug);

     ResponseVar(void);
     ~ResponseVar(void){ DBG_PRINT("ResponseVar::DTOR"); Destroy(); }
     void Destroy(void);

     void Write(FILE * pFile, int type);
     double GetInitialVal(void);
     double GetCurrentVal(void);
     UnchangeableString GetKeyword(void);
     int GetLine(void);
     int GetColumn(void);
     UnchangeableString GetFileName(void);
     UnchangeableString GetName(void);
     void SetCurrentVal(double curVal);
     void SetInitialVal(double initVal);
     char GetToken(void){ return m_Tok;}
     bool IsAugmented(void){ return m_bAug; }

     void SetLinePtr(ParameterABC * ptr){ m_pLine = ptr;}
     void SetLinePtr(TiedParamABC * ptr){ m_pTiedLine = ptr;}
     void SetColPtr(ParameterABC * ptr){ m_pCol = ptr;}
     void SetColPtr(TiedParamABC * ptr){ m_pTiedCol = ptr;}

     void WriteSim(FILE * pFile, int type);

   private :
      StringType m_Name;
      double m_InitialVal;
      double m_CurrentVal;
      StringType m_FileName;
      StringType m_Keyword;
      char m_Tok; 
      bool m_bAug;

      //Parameters from which line and column are derived
      ParameterABC * m_pLine;
      ParameterABC * m_pCol;
      //Tied parameters from which line and column are derived
      TiedParamABC * m_pTiedLine;
      TiedParamABC * m_pTiedCol;
      //Constant values of line and column
      int m_Line;
      int m_Column;
}; /* end class ResponseVar */

#endif /* RESPONSE_VAR_H */



