/******************************************************************************
File      : ResponseVarGroup.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates the response variable group, the group of response variables 
which the objective function (and possibly constraints) is based upon. Response
variables are to optimization, what observations are in regression/calibration.

Version History
05-10-04    lsm   created
01-11-05    lsm   added support for tied response variables
******************************************************************************/
#ifndef RESPONSE_VAR_GROUP_H
#define RESPONSE_VAR_GROUP_H

#include "MyHeaderInc.h"

// forward decs
class RespVarABC;
class ResponseVar;
class ValueExtractor;
class ModelBackup;
class StatsClass;
class ParameterCorrection;

/******************************************************************************
class ResponseVarGroup

  Contains a collection of response variables.
  Operations pertaining to the response variables as a whole are done here, 
  namely the operation which reads the response variables from the model output.
******************************************************************************/
class ResponseVarGroup
{
   public:
      ResponseVarGroup(void);
      ResponseVarGroup(const char * pToken);
      ~ResponseVarGroup(void){ DBG_PRINT("ResponseVarGroup::DTOR"); Destroy(); }
      void Destroy(void);

      void WriteList(FILE * pFile, int type);
      void InitializeVals(void);
      void ExtractVals(void);
      int GetNumRespVars(void);
      int GetNumTiedRespVars(void);
      RespVarABC * GetRespVarPtr(IroncladString name);
     void Write(FILE * pFile, int type);

   private:
      void InitFromFile(IroncladString respFileName);  
      void InitFromFile(IroncladString respFileName, const char * start_tag, const char * end_tag);
      void InitTiedRespVars(IroncladString pFileName);

      ResponseVar ** m_pRespVarList;
      RespVarABC ** m_pTiedRespVarList;

      //a linked list of ValueExtractor classes, one for each response file
      ValueExtractor * m_pRespFiles;

      int m_NumRespVars;
      int m_NumTiedRespVars;

   protected: //can be used by model backup class and by stats class
      RespVarABC * GetRespVarPtr(int i);
      friend class ModelBackup;
      friend class StatsClass;
      friend class ParameterCorrection;
}; /* end class ResponseVarGroup */

#endif /* RESPONSE_VAR_GROUP_H */



