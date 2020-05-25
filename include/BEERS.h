/******************************************************************************
File     : BEERS.h
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

Balanced Exploration-Exploitation Random Search (BEERS) applies a balanced
approach to randomly search a parameter space. The search initially favors
exploration of the search space and will try to maximize the distance between 
points that are evaluated. As the search progresses the algorithm favors 
exploitation and will favor points that are close to the current optimal. To
facilitate exploration vs. exploitation the algorithm maintains an archive of
every point that is ever evaluated during the search.

Version History
05-16-15    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef BEERS_H
#define BEERS_H

#include "MyHeaderInc.h"

//parent class
#include "AlgorithmABC.h"

//forward declarations
class ModelABC;

/******************************************************************************
class BEERS

******************************************************************************/
class BEERS : public AlgorithmABC
{
   public:
      BEERS(ModelABC * pModel);
      ~BEERS(void){ DBG_PRINT("BEERS::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      
      void WarmStart(void);
      int  GetCurrentIteration(void) { return m_CurSample; }

   private:
      double EstimateMaxDistance(ArchiveStruct * pA, double * pMin, double * pRange);
      void AssignModelProbs(ArchiveStruct * pA, double Fbest);      
      void CalcProbabilities(ArchiveStruct * pC, ArchiveStruct * pA, double dmax, double * pMin, double * pRange, double * pExploit, double * pExplore);
      void DestroyArchive(ArchiveStruct * pArch);
      void AdjustRanks(ArchiveStruct * pCur, ArchiveStruct * pArch);
      void WriteArchive(void);

      ModelABC * m_pModel;
      ArchiveStruct * m_pArchive;
      double * m_pMin;
      double * m_pMax;
      double * m_pRange;
      double m_MaxDist;
      double m_MinProbAccept;
      int m_NumSamples;
      int m_CurSample;
      ArchiveStruct * m_pBest;
}; /* end class BEERS */

extern "C" {
void BEERS_Program(int argC, StringType argV[]);
}

#endif /* BEERS_H */


