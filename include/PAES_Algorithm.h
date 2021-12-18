/******************************************************************************
File       : PAES_Algorithm.h
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)

PAES - Pareto archived evolution strategy 

Version History
12-29-17    lsm   created file
******************************************************************************/
#ifndef PAES_ALGORITHM_H
#define PAES_ALGORITHM_H

#include <time.h>
#include "MyTypes.h"
#include "PAES_Helpers.h"
#include "PAES_ObjFunc.h"
#include "PAES_Genetics.h"
#include "ModelABC.h"
#include "AlgorithmABC.h"

/******************************************************************************
class PAES
******************************************************************************/
class PAES : public AlgorithmABC
{
   private:
      int UpdateArchive(double * pX, int nX, double * pF, int nF);

      int    m_depth; 
      int    m_numberOfGenes;
      int    m_maximumArchiveLength;
      int    m_numberOfIterations;
      double m_mutationProbability;
      int    m_numberOfFitnessEvaluations;
      int *  m_precision;

      double m_distributionIndexForMutation;
      double m_perturbationForMutation;
  
      PAES_MutationOperator m_mutationOperator;

      FILE * m_cout;
   
      PAES_Random m_random;
      int m_seed;
      PAES_MultiobjectiveProblem * m_problem; //!< Problem to be solved  
      time_t m_startTime;
      time_t m_endTime;

      PAES_Individual * m_currentSolution;
      PAES_Individual * m_mutantSolution;
      PAES_Population * m_archiveOfSolutions;
  
      PAES_AdaptiveGrid * m_adaptiveGrid;
      int m_numberOfGridDivisions;

      int m_printFrequency;

      int m_CurIter;
      ModelABC * m_pModel;
      ArchiveStruct * m_pNonDom; //non-dominated solutions
      ArchiveStruct * m_pDom; //dominated solutions
      int m_NumNonDom; //number of non-dominated solutions
      int m_NumDom; //number of dominated solutions

      int m_pbufsize;
      int m_fbufsize;
      double * m_pbuf; //parameter buffer for parallel comm.
      double * m_fbuf; //obj. func. buffer for parallel comm.

   public:
      PAES(ModelABC * pModel);
      ~PAES(void){ DBG_PRINT("PAES::DTOR"); Destroy(); }
      void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void){ return; }
      int  GetCurrentIteration(void) { return m_CurIter; }
      void testEncodeDecode(void);
      
      // Methods
      void testAndUpdate(void);
      void initialize(PAES_MultiobjectiveProblem * problemToSolve, PAES_MutationOperator mutationOperator);
      void start(void);
      void start_parallel(int myrank, int nprocs);
      void addToArchive(PAES_Individual *solution);
      bool archiveSolution(PAES_Individual *solution);
      int  compareToArchive(PAES_Individual *solution); 
      void printStatistics(void);
      void printToFiles(char * variableFile, char * functionFile);
}; /* end class PAES */

extern "C" {
   void PAES_Program(int argC, StringType argV[]);
}


#endif /* PAES_ALGORITHM_H */

