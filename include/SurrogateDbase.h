/******************************************************************************
File     : SurrogateDbase.h
Author   : L. Shawn Matott
Copyright: 2006, L. Shawn Matott

Manages a database of model runs.

Version History
04-18-06    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef SURROGATE_DBASE_H
#define SURROGATE_DBASE_H

#include "MyHeaderInc.h"

// forward decs
class ParameterGroup;

/* defintions for the manner in which entries are added */
#define OVERWRITE_DEFAULT   (0)
#define OVERWRITE_OLDEST    (1)
#define OVERWRITE_LEAST_FIT (2)

typedef struct DBASE_ENTRY
{
   int id;
   int time_stamp;
   double run_time;
   double F;
   double * pParams;
}DbaseEntry;

/******************************************************************************
class SurrogateDbase

Manages a database of all model evaluations.
*****************************************************************************/
class SurrogateDbase
{   
   public:
      SurrogateDbase(int size, int psize, int n_models);
      ~SurrogateDbase(void){ DBG_PRINT("SurrogateDbase::DTOR"); Destroy(); }
      void Destroy(void);
      void Write(FILE * pFile);
      void Insert(ParameterGroup * pGroup, double F, int id, double run_time, int mode);
      double GetRelativeRunTime(int id);
      int GetNumStoredEvals(int id);
      int LoadBasis(int id, MyPoint * pBasis);
      void BcastEntries(ParameterGroup * pGroup, int mode);
      void BcastBestEntries(ParameterGroup * pGroup, int mode);
      double * GetBestEntry(void);
      DbaseEntry * GetBestEntry(int id);
      double GetNearestNeighbor(int id, double * pX);
      double InvDistWSSE(int id, double * pX);

   private:
      DbaseEntry * m_pDbase;
      DbaseEntry m_Temp;
      int m_CurSize;
      int m_MaxSize;
      int m_NumParams;
      int m_NumModels;
      double * m_pAvgRunTimes;
}; /* end class SurrogateDbase */

#endif /* SURROGATE_DBASE_H */


