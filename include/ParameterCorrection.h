/******************************************************************************
File     : ParameterCorrection.h
Author   : L. Shawn Matott
Copyright: 2012, L. Shawn Matott

Intereface for external parameter correction algorithm.

Version History
09-18-12    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef PARAMETER_CORRECTION_H
#define PARAMETER_CORRECTION_H

#include "MyHeaderInc.h"

// forward decs
class ParameterGroup;
class ResponseVarGroup;
class FilePair;

/******************************************************************************
class ParameterCorrection

******************************************************************************/
class ParameterCorrection
{
   public:
     ParameterCorrection(ParameterGroup * pGroup);
	 ~ParameterCorrection(void){ DBG_PRINT("ParameterCorrection::DTOR"); Destroy();}
     void Destroy(void);
          
     //misc. member functions     
     void Execute(void);
     void   WriteMetrics(FILE * pFile);
   private:
      ParameterGroup * m_pParamGroup;
      ResponseVarGroup * m_pCorrections;
      FilePair * m_FileList;
      int m_NumCorrections;
      StringType  m_ExecCmd;

      bool NearlyEqual(double a, double b);

      void SetExecCmd(IroncladString cmd);      
      void AddFilePair(FilePair * pFilePair);
}; /* end class ParameterCorrection */

#endif /* PARAMETER_CORRECTION_H */

