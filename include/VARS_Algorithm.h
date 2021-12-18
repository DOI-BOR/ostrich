/******************************************************************************
File       : VARS_Algorithm.h
Author     : L. Shawn Matott
Copyrights : 2018, L. Shawn Matott 

VARS - Variogram Analysis of Response Surfaces 

Interface class for the VARS plugin, which must be licensed separately.
 
Version History
06-18-18    lsm   created file
******************************************************************************/
#ifndef VARS_ALGORITHM_H
#define VARS_ALGORITHM_H

#include <time.h>
#include "MyTypes.h"
#include "ModelABC.h"
#include "AlgorithmABC.h"

/******************************************************************************
class VARS
******************************************************************************/
class VARS : public AlgorithmABC
{
   public:
      // required AlgorithmABC interface methods
      VARS(ModelABC * pModel);
      ~VARS(void){ DBG_PRINT("PAES::DTOR"); Destroy(); }
      void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void){ return; }
      int  GetCurrentIteration(void) { return m_CurIter; }
      
   private:
      // member functions
      int LoadPlugin(char * filename);
      int BoxCoxProgram(void);

      // member variables
      ModelABC * m_pModel;
      char * m_PluginFile;
      char * m_Subroutine;
      int m_CurIter;
}; /* end class PAES */

extern "C" {
   void VARS_Program(int argC, StringType argV[]);
}

#endif /* VARS_ALGORITHM_H */

