/******************************************************************************
File      : Observation.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates a single observation point.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
10-04-04    lsm   each observation can be assigned a different token
01-01-07    lsm   added copy CTOR and Reconfigure() routines to support Surrogate-
                  model approach
******************************************************************************/
#ifndef OBSERVATION_H
#define OBSERVATION_H

#include "MyHeaderInc.h"

#define OST_OBS_FILE "OstInterpolatedObs.txt"

/******************************************************************************
class Observation

 This class represents an observation data.
 Each observation has properties keyword, line, and column.
 The  program associates each observation with the value which is found on the 
 line and coloumn after the first occurence of the keyword in <fileName>.
******************************************************************************/
class Observation
{
   public:
      Observation(IroncladString name, double value ,double weight, 
                  IroncladString fileName, IroncladString keyword, int line,
                  int column, char tok, bool bAug, IroncladString group);

      Observation(Observation * pCopy);
      Observation(void);
      ~Observation(void){ DBG_PRINT("Observation::DTOR"); Destroy(); }
      void Destroy(void);

      void Reconfigure(IroncladString fileName, IroncladString keyword, 
                       int line, int column, char tok, bool bAug, IroncladString group);

      void Write(FILE * pFile, int type);
      void WriteSim(FILE * pFile, int type);            
      bool IsAugmented(void) { return m_bAug;}
      UnchangeableString GetKeyword(void);
      int GetLine(void);
      int GetColumn(void);     
      UnchangeableString GetFileName(void);
      UnchangeableString GetName(void);
      UnchangeableString GetGroup(void);
      void SetComputedVal(double computedVal);
      char GetToken(void){return m_Tok;}

      double CalcResidual(bool bTransformed, bool bWeighted);
      double GetMeasuredVal(bool bTransformed, bool bWeighted);
      double GetComputedVal(bool bTransformed, bool bWeighted);

   private :
      friend double GetObsWeight(Observation * pObs);

      double GetWeight(void) { return m_Weight; }

      StringType m_Name;
      double m_MeasuredVal;
      double m_ComputedVal;
      double m_Weight;
      StringType m_FileName;
      StringType m_Keyword;
      StringType m_Group;
      int m_Line;
      int m_Column;
      char m_Tok;
      bool m_bAug; //if true, include observed values in augmented OstModel file
}; /* end class Observation */

#endif /* OBSERVATION_H */



