/******************************************************************************
File     : GeneticAlgorithm.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

A Genetic Algorithm applies concepts (namely survival of the fittest and 
natural selection) from evolutionary theory to optimization problems. The 
Genetic Algorithm starts with a population of coded solutions (ChromosomePool) 
and evolves this population using the processes of Selection, Crossover and 
Mutation such that each successive genration of solutions is an improvement 
(on average) over previous generations.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
07-08-04    lsm   WriteSetup() is first action of algorithm
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
11-30-04    lsm   Replaced calls to time() with MyTime()
10-19-05    lsm   Added support for status file (used in grid computing)
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "GeneticAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "Chromosome.h"
#include "ChromosomePool.h"
#include "StatsClass.h"

#include "Exception.h"
#include "Utility.h"
#include "WriteUtility.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void GeneticAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pPopulation->SetChromosome(0, pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
GeneticAlgorithm::GeneticAlgorithm(ModelABC * pModel)
{
   RegisterAlgPtr(this);

   m_pModel = pModel;

   NEW_PRINT("ChromosomePool", 1);
   m_pPopulation = new ChromosomePool;

   MEM_CHECK(m_pPopulation);
   
   m_pStats = NULL;
   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the GA and it's member variables.
******************************************************************************/
void GeneticAlgorithm::Destroy(void)
{
   delete m_pPopulation;
   delete m_pStats;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using the GA.
******************************************************************************/
void GeneticAlgorithm::Calibrate(void)
{ 
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);

   Optimize();
   
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   if(id == 0)
   {
      sprintf(fileName, "OstOutput%d.txt", id);

      //write statistics of best parameter set to output file
      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      //write statistics of best parameter set to output file
      m_pStats->WriteStats(stdout);
   }/* end if() */
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using the GA.
******************************************************************************/
void GeneticAlgorithm::Optimize(void)
{
   StatusStruct pStatus;
   int maxGens;
   int i, id;
   double medFitness;
   //double avgFitness;
   Chromosome * pBest;
        
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   m_pPopulation->CreateComm(m_pModel);
   m_pPopulation->Initialize(); 
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }
   if(m_pModel->GetParamGroupPtr()->CheckExtraction() == true)
   {
      int num_params = m_pModel->GetParamGroupPtr()->GetNumParams();
      double * pExtract = new double[num_params];
      m_pModel->GetParamGroupPtr()->ReadParams(pExtract);
      m_pPopulation->SetChromosome(0, pExtract);
      delete [] pExtract;
   }

   maxGens = m_pPopulation->GetNumGens();
   m_StopVal = m_pPopulation->GetStopVal();  

   if(id == 0)
   {
      //write setup
      WriteSetup(m_pModel, "Real-coded Genetic Algorithm (RGA)");   
      //write banner
      WriteBanner(m_pModel, "gen    best fitness   ", " convergence value");
   }/* end if() */

   pStatus.maxIter = m_MaxGens = maxGens;
   for(i = 0; i <= maxGens; i++)
   {
      pStatus.curIter = m_CurGen = i;      
      if(IsQuit() == true){ break;}

      //evaluate current generation
      m_pPopulation->EvalFitness();

      //compute best and average
      pBest = m_pPopulation->GetBestFit();
      //avgFitness = m_pPopulation->CalcAvgFitness();
      medFitness = m_pPopulation->CalcMedianFitness();
      m_CurStop = fabs((medFitness - pBest->GetFitness())/medFitness);

      if(id == 0)
      {         
         //write results
         m_pPopulation->ConvertChromosome(pBest);
         WriteRecord(m_pModel, i, pBest->GetFitness(), m_CurStop);
         pStatus.pct  = ((float)100.00*(float)i)/(float)m_MaxGens;
         pStatus.numRuns = m_pModel->GetCounter();
         WriteStatus(&pStatus);

         //create next generation
         if(i < maxGens)
         {
            m_pPopulation->CreateNxtGen();
         }
      }/* end if() */

      if(m_CurStop < m_StopVal)
      { 
         pBest = m_pPopulation->GetBestFit();
         m_pPopulation->ConvertChromosome(pBest);
         pStatus.pct = 100.00;
         break;
      }

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */

   m_pModel->Execute();

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   {   
      WriteOptimal(m_pModel, pBest->GetFitness());
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      //write algorithm metrics
      WriteAlgMetrics(this);
   }
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out setup and metrics for the algorithm.
******************************************************************************/
void GeneticAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Real-coded Genetic Algorithm (RGA)\n");
   fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
   fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
   fprintf(pFile, "Max Generations         : %d\n", m_MaxGens);
   fprintf(pFile, "Actual Generations      : %d\n", m_CurGen);
   m_pPopulation->WriteMetrics(pFile);
   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());
   m_pModel->WriteMetrics(pFile);
   if(m_CurStop <= m_StopVal)
   {
      fprintf(pFile, "Algorithm successfully converged on a solution\n");
   }
   else
   {
      fprintf(pFile, "Algorithm failed to converge on a solution, more generations may be needed\n");
   }   
}/* end WriteMetrics() */

/******************************************************************************
GA_Program()

Calibrate the model using the GA.
******************************************************************************/
void GA_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("GeneticAlgorithm", 1);
   GeneticAlgorithm * GA = new GeneticAlgorithm(model);

   MEM_CHECK(GA);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { GA->Calibrate(); }
   else { GA->Optimize(); }

   delete GA;
   delete model;
} /* end GA_Program() */


