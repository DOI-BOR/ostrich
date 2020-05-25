/******************************************************************************
File      : ObservationGroup.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates the observation group, the group of observations which the 
objective function is based upon.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added ObsToken
08-17-04    lsm   added reporting of memory allocations, ValueExtractors are
                  now part of ObservationGroup
10-04-04    lsm   moved observation tokens to individual observations
01-01-07    lsm   Added copy CTOR and Read/Write routines to support Surrogate-
                  model approach. Added ExcludeObs() subroutine to support the
                  "hold" observations functionality.
******************************************************************************/
#ifndef OBSERVATION_GROUP_H
#define OBSERVATION_GROUP_H

#include "MyHeaderInc.h"

// forward decs
class Observation;
class ValueExtractor;

/******************************************************************************
class ObservationGroup

  contains a collection of observation points
  operations pertaining to the observation group as a whole are done here.
******************************************************************************/
class ObservationGroup
{
   public:
      ObservationGroup(void);
      //copy CTOR, used by surrogate models
      ObservationGroup(ObservationGroup * pCopy, UnmoveableString pFileName); 
      ~ObservationGroup(void){ DBG_PRINT("ObservationGroup::DTOR"); Destroy(); }
      void Destroy(void);

      void WriteList(FILE * pFile, int type);
      void ExtractVals(void);
      int GetNumObs(void);
      int GetNumGroups(void);
      UnchangeableString GetGroup(int whichGroup);
      void ReadObservations(double * obs);
      void WriteObservations(Ironclad1DArray obs);
      void Write(FILE * pFile, int type, double * F);
      Observation * GetObsPtr(IroncladString name);
      Observation * GetObsPtr(int i);
      void ExcludeObs(UnchangeableString obs);

   private:
      void InitFromFile(IroncladString obsFileName);        

      Observation ** m_pObsList;

      //a linked list of ValueExtractor classes, one for each observation file
      ValueExtractor * m_pObsFiles;

      int m_NumObs;
      int m_NumGroups;
}; /* end class ObservationGroup */

#endif /* OBSERVATION_GROUP_H */



