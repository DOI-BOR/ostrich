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

#ifndef PARAMETER_GROUP_H
#define PARAMETER_GROUP_H

#include "MyHeaderInc.h"

//forward decs
class DatabaseABC;
class ParameterABC;
class TiedParamABC;
class GeomParamABC;
class SpecialParam;
class FilePair;
class FilePipe;
struct AUG_VERT_LIST_STRUCT;
struct AUG_CIRCLE_STRUCT;

/******************************************************************************
class ParameterGroup

 Represents the collection of ingeter, continuou and combinatorial parameters 
 and deals with the group of parameters as a whole unit.
******************************************************************************/
class ParameterGroup
{
   public:
     ParameterGroup(void);
     ~ParameterGroup(void){ DBG_PRINT("ParameterGroup::DTOR"); Destroy();}
      void Destroy(void);

     void SubIntoFile(FilePipe * pPipe);
	 void SubIntoDbase(DatabaseABC * pDbase);
     void WriteDatabaseParameter(DatabaseABC * pDbase, char * find, char * replace);
     void Write(FILE * pFile, int type);     
     ParameterABC * GetParamPtr(int i);
     ParameterABC * GetParamPtr(IroncladString name);
     TiedParamABC * GetTiedParamPtr(IroncladString name);
	 //SpecialParam * GetSpecialParamPtr(IroncladString name);
     int GetNumParams(void);
     int GetNumTiedParams(void){ return m_NumTied;}
     int GetNumSpecialParams(void){ return m_NumSpecial;}
     void ReadParams(double * p);
     double WriteParams(Ironclad1DArray p);
     void CheckTemplateFiles(FilePair * pList);
     void CheckMnemonics(void);
     bool FixGeometry(void);
     void CheckBounds(void);
     void ExcludeParam(UnchangeableString prm);
     void WriteSuperMuseArgs(FILE * pFile);
	 void GetSpecialConstraints(double * pSC);
	 void ConfigureSpecialParams(double minObj, double * minCon);
	 void InitSpecialParams(IroncladString pFileName);
	 void EnableSpecialParams(void);
     bool ExtractInitialValue(char * name, bool bFixFmt, double * pVal);
     MetaParameter GetMetaParam(IroncladString name);
     bool CheckExtraction(void){ return m_bExtracted;}

   private:      
      ParameterABC ** m_pList;
      ParameterABC ** m_pExcl;
      TiedParamABC ** m_pTied;
      GeomParamABC ** m_pGeom;
	  SpecialParam ** m_pSpecial;

      int m_NumParams;
      char ** m_ParamNameList; //list of parameter names
      int m_NumTied;
      int m_NumGeom;
	  int m_NumSpecial;
      int m_NumExcl;
      void InitFromFile(IroncladString pParamFileName);
      int  CountParams(IroncladString pFileName);
      void GetParameterNames(IroncladString pFileName);
      int  GetNextEmptyParamIdx(void);
      void InitRealParams(IroncladString pFileName);
      void InitIntParams(IroncladString pFileName);
      void InitComboParams(IroncladString pFileName);
      void InitTiedParams(IroncladString pFileName);
      void InitGeomParams(IroncladString pFileName);
      bool m_bExtracted;
	  struct AUG_VERT_LIST_STRUCT * InitAugVertex(IroncladString xstr, IroncladString ystr,
                                    IroncladString zstr);
	  struct AUG_CIRCLE_STRUCT * InitAugCircle(IroncladString xstr, IroncladString ystr,
                                    IroncladString zstr, IroncladString rstr);
}; /* end class ParameterGroup */

#endif /* PARAMETER_GROUP_H */

