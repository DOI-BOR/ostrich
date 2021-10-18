/******************************************************************************
File      : ParameterGroup.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates a group of parameters. The optimization routines will attempt to
find the values of the parameters that minimize the objective function.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-19-03    lsm   Modified to support multiple stages of unit conversion.
07-05-04    lsm   added integer and combinatorial parameter support
07-08-04    lsm   added tied parameter support
11-30-04    lsm   Added support for Geometry parameters
02-03-05    lsm   Added CheckBounds()
01-01-07    lsm   Added ExcludeParam() subroutine to support the "hold" 
                  parameters functionality.
07-18-07    lsm   Added support for SuperMUSE
******************************************************************************/

#ifndef PARAMETER_GROUP_WORKER
#define PARAMETER_GROUP_WORKER

#include "MyHeaderInc.h"
#include "ParameterWorker.h"
#include "ParameterABC.h"
#include <vector>
#include <string>
#include <filesystem>
#include <iostream>

//forward decs
class DatabaseABC;
class ParameterWorker;
class FilePair;
class FilePipe;
struct AUG_VERT_LIST_STRUCT;
struct AUG_CIRCLE_STRUCT;

/******************************************************************************
class ParameterGroup

 Represents the collection of ingeter, continuou and combinatorial parameters 
 and deals with the group of parameters as a whole unit.
******************************************************************************/
class ParameterGroupWorker {
   public:
       ParameterGroupWorker();
     ~ParameterGroupWorker(void){ DBG_PRINT("ParameterGroupWorker::DTOR"); Destroy();}
      void Destroy(void);

     void SubIntoFile(FilePipe * pPipe);
	 void SubIntoDbase(DatabaseABC * pDbase);
     void WriteDatabaseParameter(DatabaseABC * pDbase, char * find, char * replace);
     void Write(FILE * pFile, int type);     
     
     ParameterWorkerABC* GetParamPtr(int i);
     ParameterWorkerABC* GetParamPtr(IroncladString name);
     int GetNumParams(void);

     void ReadParams(double * p);
     void WriteParams(Ironclad1DArray p);
     void CheckTemplateFiles(std::vector<std::vector<std::string>> m_filePairs, std::string workerDirectory);
     void CheckMnemonics(void);
     void CheckBounds(void);
     void ExcludeParam(UnchangeableString prm);
     void WriteSuperMuseArgs(FILE * pFile);

     bool ExtractValue(char * name, bool bFixFmt, double * pVal);
     MetaParameter GetMetaParam(IroncladString name);
     bool CheckExtraction(void){ return m_bExtracted;}

     
     void SetGroupValues(ParameterWorkerABC** m_pList, ParameterWorkerABC** m_pExcl, char** m_ParamNameList,
                         int numberOfParameters, int numberOfExcluded);    

     int GetNumRealParams(void) { return m_NumRealParams; };
     int GetNumIntParams(void) { return m_NumIntParams; };


   private:      
       ParameterWorkerABC** m_pList;
       ParameterWorkerABC** m_pExcl;
      

      int m_NumParams;
      int m_NumRealParams;
      int m_NumIntParams;
      char ** m_ParamNameList; //list of parameter names

      int m_NumExcl;
      void GetParameterNames(IroncladString pFileName);
      int  GetNextEmptyParamIdx(void);
      void InitRealParams(IroncladString pFileName);
      void InitIntParams(IroncladString pFileName);
      
      bool m_bExtracted;
}; 

#endif 

