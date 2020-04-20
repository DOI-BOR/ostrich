/******************************************************************************
File     : ChromosomePool.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The ChromosomePool is a container for a set of Chromosomes (coded design 
variables). The Genetic Algorithm uses ChromosomePools to store the population 
of solutions for a given generation along with the mating pool from which the 
next generation of solutions will be produced.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes
                  Added metrics collection and reporting
                  Added intermediate output (when generation is being evaluated)
11-18-04    lsm   Added convergence criteria, based on median fitness of swarm.
03-21-05    lsm   Replaced calls to rand() with MyRand()
            lsm   Added support for seeding initial population
01-01-07    lsm   Algorithms now use an abstract model base class (ModelABC).
                  Latin Hypercube Sampling has been added as an option for 
                  initializing the GA population. User's select this option by
                  including the following line in the GA section of the 
                  input file:
                     InitPopulationMethod LHS
07-16-07    lsm   Aglorithm now supports SuperMUSE
01-13-15    lsm   Added support for asynchrounous parallel
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "ChromosomePool.h"
#include "Chromosome.h"
#include "QuadTree.h"
#include "ChromosomeCommunicator.h"
#include "Gene.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "LatinHypercube.h"
#include "SuperMUSE.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"
#include "SuperMuseUtility.h"

#define APGA_DO_WORK   (101)
#define APGA_STOP_WORK (102)

/******************************************************************************
default CTOR

Assigns member variables default/NULL values and creates a prototype 
chromosome.
******************************************************************************/
ChromosomePool::ChromosomePool(void)
{
   m_pPool     = NULL;
   m_pScratch  = NULL;
   m_Fmedian = NULL;
   m_Proto = NULL;
   m_pInit = NULL;
   m_pTrees = NULL;
   m_pAssignments = NULL;
   m_TreeSize = 0;
   m_PoolSize = 50;
   m_NumInit = 0;
   m_pComm = NULL;  
   m_Generation = 0;
   m_NumSurvivors = 1;
   m_ParallelType = PARALLEL_TYPE_SYNCH;

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   m_pBuf    = NULL;

   m_pMutCount = NULL;

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()

Frees up memory associated with the ChromosomePool class.
******************************************************************************/  
void ChromosomePool::Destroy(void)
{
   int i;

   for(i = 0; i < m_PoolSize; i++) 
   {
      if(m_pPool != NULL) 
         delete m_pPool[i];
      if(m_pScratch != NULL)  
         delete m_pScratch[i];
   }
   delete [] m_pPool;
   delete [] m_pTrees;

   for(i = 0; i < m_NumInit; i++) 
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;

   delete m_pComm;
   delete m_Proto;
   delete [] m_pBigBuf;
   delete [] m_pTmpBuf;
   delete [] m_pMyBuf; 
   delete [] m_pBuf;  
   delete [] m_pMutCount;
   delete [] m_Fmedian;
   delete [] m_pAssignments;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
TourneySelection()

Determines the mating pool by randomly selecting two chromosomes and comparing
their fitness values. The chromosome with the better fitness gains the right 
to pass its genes into the next generation. The configuration variable 
m_NumSurvivors is used to guarantee that the top chromosomes survive unchanged 
into the next generation.

The input argument specifies the number of combatants in the tournament. For
'standard' GA this is set equal to 2. For computation constrained, the number
of combatants increases as number of generations increases.
******************************************************************************/  
void ChromosomePool::TourneySelection(int nCombatants)
{
   int i, j; 
   int r1;
   int r2;
   Chromosome * pPlay1; //player one of tournament
   Chromosome * pPlay2; //player two of tournament
   Chromosome * pMax;   //chromosome with max. fitness
   double fit1;
   double fit2;
   double maxFit;
   double lastMax;

   /*-------------------------------------------
   Reserve the top m_NumSurvivors chromosomes. 
   --------------------------------------------*/
   lastMax = NEARLY_HUGE;
   pMax = NULL;
   for(i = 0; i < m_NumSurvivors; i++)
   {
      maxFit = -NEARLY_HUGE;

      for(j = 0; j < m_PoolSize; j++)
      {
         pPlay1 = m_pPool[j];
         fit1 = pPlay1->GetFitness();

         if((fit1 > maxFit) && (fit1 <= lastMax) && (pMax != pPlay1))
         {
            pMax = pPlay1;
            maxFit = fit1;
         }/* end if() */
      }/* end for() */      

      //propagate nth max. to next generation
      lastMax = maxFit;
      m_pScratch[i]->Copy(pMax);
      //printf("Max. #%d = %E\n", i, lastMax);
   }/* end for() */

   /*-------------------------------------------
   Use n-member tourney to select the remaining 
   chromosomes.
   --------------------------------------------*/
   for(i = m_NumSurvivors; i < m_PoolSize;i++)
   {
      // pick random chromosomes 
      r1 = MyRand() % m_PoolSize;
      pPlay1 = m_pPool[r1];
      fit1 = pPlay1->GetFitness();

      for(j = 0; j < (nCombatants - 1); j++)
      {
         r2 = MyRand() % m_PoolSize;
         pPlay2 = m_pPool[r2];

         //evaluate their fitness
         fit2 = pPlay2->GetFitness();

         //the better one gets to go to the nextGeneration
         if (fit2>fit1)
         {
            pPlay1 = pPlay2;
            fit1 = fit2;
         }
      }
   
      m_pScratch[i]->Copy(pPlay1); 
   } /* end for() */

   /*------------------------------
   copy the latest generation from
   scratch into the chromsome pool.
   -------------------------------*/
   for(i = 0; i < m_PoolSize; i++) 
   {       
      m_pPool[i]->Copy(m_pScratch[i]);
   }/* end for() */
} /* end TourneySelection() */

/******************************************************************************
Crossover()

Crosses over each chromsome of the population with the next one in the 
population, except those that are in the top <m_NumSurvivors>.
******************************************************************************/  
void ChromosomePool::Crossover(void)
{   
   int i;
   Chromosome * pFirst;
   Chromosome * pMom;
   Chromosome * pPop;

   //save first for later
   pMom   = m_pPool[m_NumSurvivors];
   m_pScratch[0]->Copy(pMom);
   pFirst = m_pScratch[0];

   //crossover everyone with their neighbor
   for(i = m_NumSurvivors; i < (m_PoolSize - 1); i++)
   {
      pMom = m_pPool[i];
      pPop = m_pPool[i+1];
      pMom->Crossover(pPop);  
   }/* end for() */

   //crossover last and first
   pMom = m_pPool[m_PoolSize - 1];
   pMom->Crossover(pFirst);
}/* end Crossover() */

/******************************************************************************
Mutate()

Mutates individual chromsomes of the population according to a  
pre-established mutation rate.
******************************************************************************/  
void ChromosomePool::Mutate(void)
{        
   int i;

   for(i = m_NumSurvivors; i < m_PoolSize; i++) 
   { 
      m_pPool[i]->Mutate(m_pMutCount); 
   }
}/* end Mutate() */

/******************************************************************************
CreateNxtGen()

Creates the next generation of the chromosome population using tourney selection,
crossover, and mutation.
******************************************************************************/  
void  ChromosomePool::CreateNxtGen(void)
{
  m_Generation++;
  TourneySelection(2);
  Crossover();
  Mutate();
}/* end CreateNxtGen() */

/******************************************************************************
CreateNxtGen()

Creates the next generation of the chromosome population using tourney selection,
crossover, and mutation. Adapts the GA operators as optimization proceeds.
******************************************************************************/
void  ChromosomePool::CreateNxtGen(double pct)
{
  int np = m_pComm->GetParamGroupPtr()->GetNumParams();
  int ng = m_NumGenerations;
  int nCombatants = (int)(0.5+(2.00+pct*0.5*(ng-2.00)));
  m_Generation++;
  TourneySelection(nCombatants);
  Crossover();  
  //adjust mutation rate
  for(int i = 0; i < m_PoolSize; i++) 
  { 
     m_pPool[i]->SetMutationRate(0.15*(1.00-pct)); 
  }
  Mutate();
  //freeze a certain number of genes at their optimal values
  //more are more are frozen as the optimization proceeds
  //int nFreeze = (int)(pct*(double)np);
  //FreezeGenes(nFreeze);
}/* end CreateNxtGen() */

/******************************************************************************
FreezeGenes()

For each chromosome, randomly freeze the given number of genes so that they 
are at the current global optimal.
******************************************************************************/
void ChromosomePool::FreezeGenes(int numFreeze)
{
   int np = m_pComm->GetParamGroupPtr()->GetNumParams();

   for(int i = m_NumSurvivors; i < m_PoolSize; i++)
   {
      SampleWithReplacement(-2, np);
      for(int j = 0; j < numFreeze; j++)
      {
         int k = SampleWithReplacement(1, np);
         m_pPool[i]->GetGenePtr(k)->SetValue(GetBestFit()->GetGenePtr(k)->GetValue());
      }
   }
}/* FreezeGenes() */

/******************************************************************************
CalcAvgFitness()

Calculates and returns the average fitness of the population. This parameter 
is used in the termination criteria of the genetic algorithm.
******************************************************************************/  
double ChromosomePool::CalcAvgFitness(void)
{
   double sum, avg;
   int i;
    
   sum = 0.00;

   for(i = 0; i < m_PoolSize; i++) { sum += m_pPool[i]->GetFitness(); }  
  
   avg = sum / (double)m_PoolSize;

   return avg;
} /* end CalcAvgFitness() */

/******************************************************************************
CalcMedianFitness()

Calculates and returns the median fitness of the population. This parameter 
is used in the termination criteria of the genetic algorithm.
******************************************************************************/  
double ChromosomePool::CalcMedianFitness(void)
{   
   double med;
   int i;
    
   for(i = 0; i < m_PoolSize; i++) { m_Fmedian[i] = m_pPool[i]->GetFitness(); }  
   med = CalcMedian(m_Fmedian, m_PoolSize);

   return med;
} /* end CalcMedianFitness() */

/******************************************************************************
GetBestFit()

Retrieves the chromosome that has the best fitness value.
******************************************************************************/
Chromosome * ChromosomePool::GetBestFit(void)
{
   double bestFitVal;
   double tmp;
   int bestFitIdx;
   int i;
   
   bestFitIdx = 0;
   bestFitVal = m_pPool[0]->GetFitness();   
   
   for(i = 0; i < m_PoolSize; i++)
   {
      tmp = m_pPool[i]->GetFitness();  
      if(tmp > bestFitVal)
      {
         bestFitVal = tmp;
         bestFitIdx = i;
      }/* end if() */
   }/* end for() */
    
   return (m_pPool[bestFitIdx]);   
}/* end GetBestFit() */

/******************************************************************************
EvalFitness()

Evaluates the fitness of each chromosome in the pool.
******************************************************************************/
void ChromosomePool::EvalFitness(void)
{   
   int i, n, id;   

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   
   if(n == 1) //serial or SuperMUSE-parallel execution
   {
      if(IsSuperMUSE() == false)
      {
         WriteInnerEval(WRITE_GA, m_PoolSize, '.');
         for(i = 0; i < m_PoolSize; i++) 
         { 
            WriteInnerEval(i+1, m_PoolSize, '.');
            m_pComm->MakeParameterCorrections(m_pPool[i]); 
            m_pComm->EvalFitness(m_pPool[i]); 
         }
         WriteInnerEval(WRITE_ENDED, m_PoolSize, '.');
      }
      else //SuperMUSE
      {
         EvalFitSuperMUSE();
      }
   }/* end if() */
   else /* MPI-parallel execution */
   {
      if(m_ParallelType == PARALLEL_TYPE_SYNCH)
      {
         if(id == 0)
         {
            for(i = 0; i < m_PoolSize; i++) 
            { 
               m_pComm->MakeParameterCorrections(m_pPool[i]); 
            }
         }
         BcastPopulation();
         EvalFitParallel();
      }/* end if() */
      else //(m_ParallelType == PARALLEL_TYPE_ASYNCH)
      {
         EvalFitnessAsynch(id, n);
      }/* end else() */
   }/* end else() */
} /* end EvalFitness() */

/******************************************************************************
EvalFitnessAsynch()

When in asynchronous parallel, master sends each parameter set out to first
available slave.
******************************************************************************/
void ChromosomePool::EvalFitnessAsynch(int rank, int nprocs)
{
   static double fbest = HUGE_VAL;
   int i, j;
   MPI_Status mpi_status;
   int ii, signal, num, sid, num_recv, nstops;
   double f;
   bool bDone = false;
   ParameterGroup * pGroup =  m_pComm->GetParamGroupPtr();

   num = pGroup->GetNumParams();

   //allocate space for data message
   if(m_pMyBuf == NULL)
   {   
      NEW_PRINT("double", num);
      m_pMyBuf = new double[num];
      MEM_CHECK(m_pMyBuf);
   }

   //allocate space for slave assignments
   if(m_pAssignments == NULL)
   {   
      NEW_PRINT("int", nprocs);
      m_pAssignments = new int[nprocs];
      MEM_CHECK(m_pAssignments);
   }

   nstops = 0;

   if(rank == 0)
   {
      //adjust parameter values using rules engine
      for(i = 0; i < m_PoolSize; i++) 
      { 
         m_pComm->MakeParameterCorrections(m_pPool[i]); 
      }

      /*------------------------------------------------
      Send initial parameter sets off to waiting slaves
      -------------------------------------------------*/

      //output banner
      WriteInnerEval(WRITE_GA, m_PoolSize, '.');

      //assign initial work to slaves
      for(i = 1; i < nprocs; i++)
      {
         if(i <= m_PoolSize)
         {
            m_pAssignments[i] = i-1;

            //fill data message
            for(j = 0; j < num; j++)
            {
               m_pMyBuf[j] = m_pPool[i-1]->GetGenePtr(j)->GetValue();
            }/* end for() */

            // send work to slave
            signal = APGA_DO_WORK;
            MPI_Send(&signal,1,MPI_INT,i,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
            MPI_Send(&(m_pMyBuf[0]), num, MPI_DOUBLE, i, MPI_DATA_TAG, MPI_COMM_WORLD);
         }
         else
         {
            signal = APGA_STOP_WORK;
            MPI_Send(&signal,1,MPI_INT,i,MPI_REQUEST_TAG,MPI_COMM_WORLD);
            nstops++;
         }
      }/* end for(each slave) */
   }/* end if() */

   num_recv = 0;
   while(bDone == false)
   {
      if(rank == 0)
      {
         //receive result from slave and process
         MPI_Recv(&f, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_RESULTS_TAG, MPI_COMM_WORLD, &mpi_status);
         num_recv++;
         sid = mpi_status.MPI_SOURCE;
         WriteInnerEval(num_recv, m_PoolSize, '.');
         ii = m_pAssignments[sid];
         m_pPool[ii]->SetFitness(-f);

         if(f < fbest)
         {
            fbest = f;
            SaveModel(sid);
         }/* end if() */

         //assign more work
         if(i <= m_PoolSize)
         {
            m_pAssignments[sid] = i-1;

            //fill data message
            for(j = 0; j < num; j++)
            {
               m_pMyBuf[j] = m_pPool[i-1]->GetGenePtr(j)->GetValue();
            }/* end for() */

            // send work to slave
            signal = APGA_DO_WORK;
            MPI_Send(&signal,1,MPI_INT,sid,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
            MPI_Send(&(m_pMyBuf[0]), num, MPI_DOUBLE, sid, MPI_DATA_TAG, MPI_COMM_WORLD);
            i++;
         }
         else // send stop work message to the slave
         {
            signal = APGA_STOP_WORK;
            MPI_Send(&signal,1,MPI_INT,sid,MPI_REQUEST_TAG,MPI_COMM_WORLD);
            nstops++;
            if(nstops == (nprocs-1)) //quit when all slaves have been told to stop
            {
               WriteInnerEval(WRITE_ENDED, m_PoolSize, '.');
               bDone = true;
            }
         }
      }/* end if(master) */
      else /* slave processing */
      {
         MPI_Recv(&signal,1,MPI_INT,0,MPI_REQUEST_TAG,MPI_COMM_WORLD, &mpi_status); 
         if(signal == APGA_DO_WORK)
         {
            num_recv++;
            MPI_Recv(&(m_pMyBuf[0]), num, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
            pGroup->WriteParams(m_pMyBuf);            
            f = RunModel();
            MPI_Send(&f, 1, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);
         }/* end if() */
         else
         {
            bDone = true;
         }
      }/* end else(slave) */
   }/* end while() */

   //synch up processors
   MPI_Barrier(MPI_COMM_WORLD);
}/* end EvalFitnessAsynch() */

/******************************************************************************
BcastPopulation()

When in parallel, only the master computes the random processes of tourney 
selection, crossover and mutation. All the other processors just compute the
fitness functions. The BcastPopulation() routine is called upon to broadcast 
the current population members from the master processor to all of the slave 
processors.
******************************************************************************/
void ChromosomePool::BcastPopulation(void)
{   
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix
   pop_size = m_PoolSize;
   num_vars = m_pPool[0]->GetNumGenes();

   buf_size = pop_size*num_vars;
   if(m_pBuf == NULL)
   {
      NEW_PRINT("double", buf_size);
      m_pBuf = new double[buf_size];   
      MEM_CHECK(m_pBuf);
   }

   for(i = 0; i < buf_size; i++){ m_pBuf[i] = 999.99;}

   //fill up the flattened matrix
   for(i = 0; i < num_vars; i++)
   {
      for(j = 0; j < pop_size; j++)
      {
         idx = (num_vars)*j + i;
         m_pBuf[idx] = m_pPool[j]->GetGenePtr(i)->GetValue();
      }/* end for() */
   }/* end for() */

   //broadcast the flattened matrix
   MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //use the flattened matrix to fill gene pool
   for(i = 0; i < num_vars; i++)
   {
      for(j = 0; j < pop_size; j++)
      {
         idx = (num_vars)*j + i;
         m_pPool[j]->GetGenePtr(i)->SetValue(m_pBuf[idx]);
      }/* end for() */
   }/* end for() */
}/* end BcastPopulation() */

/******************************************************************************
EvalFitParallel()

Compute fitness of entire population in parallel using MPI. Each processor 
evaluates a predetermined number of population members, based on their 
processor id.
******************************************************************************/
void ChromosomePool::EvalFitParallel(void)
{
   int i ,j, num_procs, id, bufsize, idx;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   //allocate space for intermediate buffers, if necessary
   bufsize = (m_PoolSize/num_procs) + 1;
   if(m_pMyBuf == NULL)
   {   
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];
      NEW_PRINT("double", bufsize);
      m_pTmpBuf = new double[bufsize];
      NEW_PRINT("double", m_PoolSize);
      m_pBigBuf = new double[m_PoolSize];   
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   for(i = 0; i < m_PoolSize; i++) 
   { 
      if((i % num_procs) == id)
      { 
         m_pComm->EvalFitness(m_pPool[i]);         
         m_pMyBuf[j] = m_pPool[i]->GetFitness();
         m_pTmpBuf[j] = m_pMyBuf[j];
         j++;
      }/* end if() */
   }/* end for() */

   //gather results
   for(i = 0; i < num_procs; i++)
   {
      //receive someones buf, this will clobber myBuf
      MPI_Bcast(m_pMyBuf, bufsize, MPI_DOUBLE, i, MPI_COMM_WORLD);

      for(j = 0; j < bufsize; j++)
      {
         idx = (num_procs * j) + i; //idx maps myBuf into bigBuf

         if(idx < m_PoolSize)
         {
            m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            m_pMyBuf[j] = m_pTmpBuf[j]; //restore myBuf...clobbered by bcast
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //stuff results into population
   for(i = 0; i < m_PoolSize; i++)
   {
      m_pPool[i]->SetFitness(m_pBigBuf[i]);
   }/* end for() */
}/* end EvalFitParallel() */

/******************************************************************************
EvalFitSuperMUSE()

Compute fitness of entire population using SuperMUSE. This routine interfaces
with the RepeatTasker SuperMUSE program, which assigns model evaluations to
SuperMUSE clients on a first-come-first-served basis.
******************************************************************************/
void ChromosomePool::EvalFitSuperMUSE(void)
{  
   double val;
   bool bOk; 
   int pop_size;   
   int i;
   ParameterGroup * pGroup;
   SuperMUSE * pSMUSE = GetSuperMusePtr();

   pop_size = m_PoolSize;

   /* ----------------------------------------------------------------
   Generate task file that describes the desired parallel evaluations.
   This is somewhat analogous to the BcastPopulation() operation used
   for MPI-parallel operations. Write the parameter values of each 
   population member as entries in the task file.

   Entries are first accumlated into a temp file to prevent the 
   SuperMUSE RepeatTasker program from prematurely processing the task
   file.
   ---------------------------------------------------------------- */
   for(i = 0; i < pop_size; i++)
   {
      //stuff the parameter group with values
      pGroup = ConvertChromosome(m_pPool[i]);
         
      //pass group to supermuse
      pSMUSE->WriteTask(pGroup);
   }/* end for() */

   //Finish task file (this will cause RepeatTasker to begin processing the job)
   pSMUSE->FinishTaskFile();

   //wait for SuperMUSE to report back (via the success or error files)
   bOk = pSMUSE->WaitForTasker();

   if(bOk == false) //SuperMUSE failed
   {
      LogError(ERR_SMUSE, "Reverting to serial execution.");
      DisableSuperMUSE();
      EvalFitness();
   }
   else //SuperMUSE was successful
   {
      for(i = 0; i < pop_size; i++)
      {
         /* -----------------------------------------------
         Stuff the parameter group with ith population
         member. This ensures that each objective function 
         gets associated with the correct parameter values.
         ------------------------------------------------ */
         pGroup = ConvertChromosome(m_pPool[i]);

         //stuff i-th result into chromosome pool
         val = pSMUSE->GatherResult(i);
         m_pPool[i]->SetFitness(-val);
      }/* end for() */
   }/* end else() */
}/* end EvalFitSuperMUSE() */

/******************************************************************************
Initialize()

Initializes the population. First, all parameters are assigned default values
and then the user input file is checked for overriding values.
******************************************************************************/
void ChromosomePool::Initialize(void)
{   
   FILE * pFile;
   int popSize, num;
   double rate, lwr, upr;
   double * pVals;
   char * pTok;
   LatinHypercube * pLHS = NULL;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   int i, j, k, lvl, idx;
   IroncladString name = GetInFileName();
   
   m_Proto = NULL;

   //read in population size and mutation rate   
   m_NumGenerations = 10;
   m_NumInit = 0;
   popSize = 50;
   rate = 0.05;
   m_NumSurvivors = 1;
   m_InitType = RANDOM_INIT;
   m_StopVal = 0.0001;

   pFile = fopen(name, "r");
   if(pFile != NULL) 
   {
      if(CheckToken(pFile, "BeginGeneticAlg", name) == true)
      {
         FindToken(pFile, "EndGeneticAlg", name);
         rewind(pFile);

         FindToken(pFile, "BeginGeneticAlg", name);
         line = GetNxtDataLine(pFile, name);
      
         while(strstr(line, "EndGeneticAlg") == NULL)
         {
            if(strstr(line, "ParallelMethod") != NULL)
            {
               sscanf(line, "%s %s", tmp, tmp2);
               MyStrLwr(tmp2);
               if(strcmp(tmp2, "synchronous") == 0)
               {
                  m_ParallelType = PARALLEL_TYPE_SYNCH;
               }
               else if(strcmp(tmp2, "asynchronous") == 0)
               {
                  m_ParallelType = PARALLEL_TYPE_ASYNCH;
               }
               else
               {
                  m_ParallelType = PARALLEL_TYPE_SYNCH;
               }
            }
            if(strstr(line, "PopulationSize") != NULL)
            {
               sscanf(line, "%s %d", tmp, &popSize);
            }
            else if(strstr(line, "MutationRate") != NULL)
            {
               sscanf(line, "%s %lf", tmp, &rate);
            }
            else if(strstr(line, "Survivors") != NULL)
            {
               sscanf(line, "%s %d", tmp, &m_NumSurvivors);
            }
            else if(strstr(line, "NumGenerations") != NULL)
            {
               sscanf(line, "%s %d", tmp, &m_NumGenerations);
            }
            else if(strstr(line, "InitPopulationMethod") != NULL)
            {
               sscanf(line, "%s %s", tmp, tmp2);
               MyStrLwr(tmp2);
               if(strcmp(tmp2, "random") == 0) {m_InitType = RANDOM_INIT;}
               else if(strcmp(tmp2, "quadtree") == 0) {m_InitType = QUAD_TREE_INIT;}
               else if(strcmp(tmp2, "lhs") == 0) {m_InitType = LHS_INIT;}
            }
            else if(strstr(line, "ConvergenceVal") != NULL)
            {
               sscanf(line, "%s %lf", tmp, &m_StopVal);
            }
            line = GetNxtDataLine(pFile, name);
         }/* end while() */
      }/* end if() */
      else
      {
         LogError(ERR_FILE_IO, "Using default algorithm setup.");
      }/* end else() */

      /* initialize some or all pop. members to specied values */
      rewind(pFile);
      if(CheckToken(pFile, "BeginInitParams", name) == true)
      {
         FindToken(pFile, "EndInitParams", name);
         rewind(pFile);

         //allocate space for the parameter list
         num = m_pComm->GetParamGroupPtr()->GetNumParams();

         //count the number of entries
         FindToken(pFile, "BeginInitParams", name);
         line = GetNxtDataLine(pFile, name);
         m_NumInit = 0;
         while(strstr(line, "EndInitParams") == NULL)
         {
            m_NumInit++;
            line = GetNxtDataLine(pFile, name);
         }/* end while() */

         //allocate space for entries
         if(m_NumInit > 0)
         {
            NEW_PRINT("double *", m_NumInit);
            m_pInit = new double * [m_NumInit];
            MEM_CHECK(m_pInit);
            for(i = 0; i < m_NumInit; i++)
            { 
               NEW_PRINT("double", num);
               m_pInit[i] = new double[num];
               MEM_CHECK(m_pInit[i]);
            }
         }/* end if() */

         //read in entries
         rewind(pFile);
         FindToken(pFile, "BeginInitParams", name);
         line = GetNxtDataLine(pFile, name);
         i = 0;
         while(strstr(line, "EndInitParams") == NULL)
         {
            pTok = line;
            //extract values, one-by-one, making any necessary conversions
            for(k = 0; k < num; k++)
            {
               j = ExtractString(pTok, tmp);
               j = ValidateExtraction(j, k, num, "ChromosomePool::Initialize()");
               pTok += j;            
               m_pInit[i][k] = m_pComm->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
            }/* end for() */                  
            i++;
            line = GetNxtDataLine(pFile, name);
         }/* end while() */
      }/* end if() */

      fclose(pFile);
   }/* end if() */   
   
   //check population and mutation rate
   if(popSize <= 0)
   {
      LogError(ERR_FILE_IO, "Invalid population size");
      ExitProgram(1);
   }/* end if() */
   if((rate < 0.00) || (rate > 1.00))
   {
      LogError(ERR_FILE_IO, "Invalid mutation rate");
      ExitProgram(1);
   }/* end if() */
   if(m_NumGenerations <= 0)
   {
      LogError(ERR_FILE_IO, "Invalid number of generations");
      ExitProgram(1);
   }
   
   m_Generation = 0;
   m_Proto = m_pComm->CreateProto(rate);

   //init. mutation counters
   NEW_PRINT("int", m_Proto->GetNumGenes());
   m_pMutCount = new int[m_Proto->GetNumGenes()];
   MEM_CHECK(m_pMutCount);
   for(i = 0; i < m_Proto->GetNumGenes(); i++){ m_pMutCount[i] = 0;}
   
   m_PoolSize = popSize;
   NEW_PRINT("Chromosome *", m_PoolSize);
   m_pPool = new Chromosome * [m_PoolSize];   
   MEM_CHECK(m_pPool);

   NEW_PRINT("Chromosome *", m_PoolSize);
   m_pScratch = new Chromosome * [m_PoolSize];   
   MEM_CHECK(m_pScratch);

   NEW_PRINT("double", m_PoolSize);
   m_Fmedian = new double [m_PoolSize];   
   MEM_CHECK(m_Fmedian);

   if(m_InitType == LHS_INIT)
   {
      NEW_PRINT("double", m_Proto->GetNumGenes());
      pVals = new double[m_Proto->GetNumGenes()];
      MEM_CHECK(pVals);

      NEW_PRINT("LatinHypercube", 1);
      pLHS = new LatinHypercube(m_Proto->GetNumGenes(), m_PoolSize);
      MEM_CHECK(pLHS);

      for(j = 0; j < m_Proto->GetNumGenes(); j++)
      { 
          lwr = m_Proto->GetGenePtr(j)->GetLwr();
          upr = m_Proto->GetGenePtr(j)->GetUpr();
          pLHS->InitRow(j, lwr, upr);
      }/* end for() */
   }/* end if() */

   lvl = idx = 0;
   for(i = 0; i < m_PoolSize; i++)
   {
      if(m_InitType == RANDOM_INIT)
      {
         m_pPool[i]    = m_Proto->CreateRandomChromo();
         m_pScratch[i] = m_Proto->CreateRandomChromo();
      } /* end if() */
      else if(m_InitType == QUAD_TREE_INIT)
      {
         //initialize quad trees if needed
         if(m_pTrees == NULL)
         {
            m_TreeSize = m_Proto->GetNumGenes();
            NEW_PRINT("QuadTree", m_TreeSize);
            m_pTrees = new QuadTree[m_TreeSize];            
            for(j = 0; j < m_TreeSize; j++)
            { 
               lwr = m_Proto->GetGenePtr(j)->GetLwr();
               upr = m_Proto->GetGenePtr(j)->GetUpr();
               m_pTrees[j].Init(lwr, upr);
            }/* end for() */
         }/* end if() */

         pVals = GetTreeCombo(lvl, idx, m_pTrees, m_TreeSize);
         //expand tree if needed.
         if(pVals == NULL)
         {
            for(j = 0; j < m_TreeSize; j++){ m_pTrees[j].Expand();}
            lvl++;
            idx = 0;            
            pVals = GetTreeCombo(lvl, idx, m_pTrees, m_TreeSize);
         }
         idx++;
         m_pPool[i]    = m_Proto->CreateChromo(pVals);
         m_pScratch[i] = m_Proto->CreateChromo(pVals);
         delete [] pVals;
      }/* end else if(QUAD_TREE_INIT) */
      else if(m_InitType == LHS_INIT)
      {
         for(j = 0; j < m_Proto->GetNumGenes(); j++){ pVals[j] = pLHS->SampleRow(j);}
         m_pPool[i]    = m_Proto->CreateChromo(pVals);
         m_pScratch[i] = m_Proto->CreateChromo(pVals);
      }/* end else() */
   }/* end for() */

   //seed initial population
   for(i = 0; i < m_NumInit; i++)
   {
      delete m_pPool[i];
      m_pPool[i] = m_Proto->CreateChromo(m_pInit[i]);
      delete m_pScratch[i];
      m_pScratch[i] = m_Proto->CreateChromo(m_pInit[i]);
   }

   if(m_InitType == LHS_INIT)
   {
      delete pLHS;
      delete [] pVals;
   }

   ((ModelChromoComm *)(m_pComm))->SetMaxEvals((m_NumGenerations+1)*m_PoolSize);
} /* end Initialize() */

/******************************************************************************
SetChromosome()

Replace the ith chromosome with the given vector.
******************************************************************************/
void ChromosomePool::SetChromosome(int i, double * vals)
{
   delete m_pPool[i];
   m_pPool[i] = m_Proto->CreateChromo(vals);
   delete m_pScratch[i];
   m_pScratch[i] = m_Proto->CreateChromo(vals);   
}/* end SetChromosome() */

/******************************************************************************
Initialize()

Initializes the population so that it is on a budget.
******************************************************************************/
void ChromosomePool::Initialize(int * budget)
{   
   FILE * pFile;
   int popSize, num, np;
   double rate, lwr, upr;
   double * pVals;
   char * pTok;
   LatinHypercube * pLHS = NULL;
   char * line;
   char tmp[DEF_STR_SZ];
   int i, j, k, lvl, idx;
   IroncladString name = GetInFileName();
   
   //assign default GA parameters
   np = m_pComm->GetParamGroupPtr()->GetNumParams();
   *budget = 1000;
   m_NumGenerations = (int)(0.5 + 2*np + sqrt((double)*budget));
   m_NumInit = 0;
   popSize = (int)(0.5 + (*budget/m_NumGenerations));
   rate = 0.15; //initial mutation rate is high (15%)
   m_NumSurvivors = (int)(MyMax(1.00, 0.5+0.05*popSize));
   m_InitType = LHS_INIT;
   m_StopVal = -1.00;

   pFile = fopen(name, "r");
   if(pFile != NULL) 
   {
      if(CheckToken(pFile, "BeginGeneticAlg", name) == true)
      {
         FindToken(pFile, "EndGeneticAlg", name);
         rewind(pFile);

         FindToken(pFile, "BeginGeneticAlg", name);
         line = GetNxtDataLine(pFile, name);
      
         while(strstr(line, "EndGeneticAlg") == NULL)
         {
            if(strstr(line, "Budget") != NULL)
            {
               sscanf(line, "%s %d", tmp, budget);
               if(*budget <= 0) *budget = 1000;
            }
            line = GetNxtDataLine(pFile, name);
         }/* end while() */
      }/* end if() */
      else
      {
         LogError(ERR_FILE_IO, "Using default algorithm setup.");
      }/* end else() */

      /* initialize some or all pop. members to specied values */
      rewind(pFile);
      if(CheckToken(pFile, "BeginInitParams", name) == true)
      {
         FindToken(pFile, "EndInitParams", name);
         rewind(pFile);

         //allocate space for the parameter list
         num = m_pComm->GetParamGroupPtr()->GetNumParams();

         //count the number of entries
         FindToken(pFile, "BeginInitParams", name);
         line = GetNxtDataLine(pFile, name);
         m_NumInit = 0;
         while(strstr(line, "EndInitParams") == NULL)
         {
            m_NumInit++;
            line = GetNxtDataLine(pFile, name);
         }/* end while() */

         //allocate space for entries
         if(m_NumInit > 0)
         {
            NEW_PRINT("double *", m_NumInit);
            m_pInit = new double * [m_NumInit];
            MEM_CHECK(m_pInit);
            for(i = 0; i < m_NumInit; i++)
            { 
               NEW_PRINT("double", num);
               m_pInit[i] = new double[num];
               MEM_CHECK(m_pInit[i]);
            }
         }/* end if() */

         //read in entries
         rewind(pFile);
         FindToken(pFile, "BeginInitParams", name);
         line = GetNxtDataLine(pFile, name);
         i = 0;
         while(strstr(line, "EndInitParams") == NULL)
         {
            pTok = line;
            //extract values, one-by-one, making any necessary conversions
            for(k = 0; k < num; k++)
            {
               j = ExtractString(pTok, tmp);
               j = ValidateExtraction(j, k, num, "ChromosomePool::Initialize()");
               pTok += j;            
               m_pInit[i][k] = m_pComm->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
            }/* end for() */                  
            i++;
            line = GetNxtDataLine(pFile, name);
         }/* end while() */
      }/* end if() */

      fclose(pFile);
   }/* end if() */   

   /* -------------------------------------------------------------------
   Adjust population size and max. gens to reflect user-defined budget
   ------------------------------------------------------------------- */
   if(*budget > popSize*m_NumGenerations)
   {
      //inc. max. gens
      m_NumGenerations = *budget/popSize;
   }
   else if(*budget < (popSize*m_NumGenerations))
   {
      if(*budget < (popSize*3)) //revise pop. size
      {
         popSize = (int)MyMax(*budget/3, 3);
         m_NumGenerations = *budget/popSize;
      }
      else //revise max gens
      {
         m_NumGenerations = *budget/popSize;
      }
   }
   if(popSize*m_NumGenerations < *budget) m_NumGenerations++;
      
   m_Generation = 0;
   m_Proto = m_pComm->CreateProto(rate);

   //init. mutation counters
   NEW_PRINT("int", m_Proto->GetNumGenes());
   m_pMutCount = new int[m_Proto->GetNumGenes()];
   MEM_CHECK(m_pMutCount);
   for(i = 0; i < m_Proto->GetNumGenes(); i++){ m_pMutCount[i] = 0;}
   
   m_PoolSize = popSize;
   NEW_PRINT("Chromosome *", m_PoolSize);
   m_pPool = new Chromosome * [m_PoolSize];   
   MEM_CHECK(m_pPool);

   NEW_PRINT("Chromosome *", m_PoolSize);
   m_pScratch = new Chromosome * [m_PoolSize];   
   MEM_CHECK(m_pScratch);

   NEW_PRINT("double", m_PoolSize);
   m_Fmedian = new double [m_PoolSize];   
   MEM_CHECK(m_Fmedian);

   NEW_PRINT("double", m_Proto->GetNumGenes());
   pVals = new double[m_Proto->GetNumGenes()];
   MEM_CHECK(pVals);

   NEW_PRINT("LatinHypercube", 1);
   pLHS = new LatinHypercube(m_Proto->GetNumGenes(), m_PoolSize);
   MEM_CHECK(pLHS);

   for(j = 0; j < m_Proto->GetNumGenes(); j++)
   { 
         lwr = m_Proto->GetGenePtr(j)->GetLwr();
         upr = m_Proto->GetGenePtr(j)->GetUpr();
         pLHS->InitRow(j, lwr, upr);
   }/* end for() */

   lvl = idx = 0;
   for(i = 0; i < m_PoolSize; i++)
   {
      for(j = 0; j < m_Proto->GetNumGenes(); j++){ pVals[j] = pLHS->SampleRow(j);}
      m_pPool[i]    = m_Proto->CreateChromo(pVals);
      m_pScratch[i] = m_Proto->CreateChromo(pVals);
   }/* end for() */

   //seed initial population
   for(i = 0; i < m_NumInit; i++)
   {
      delete m_pPool[i];
      m_pPool[i] = m_Proto->CreateChromo(m_pInit[i]);
      delete m_pScratch[i];
      m_pScratch[i] = m_Proto->CreateChromo(m_pInit[i]);
   }
   delete pLHS;
   delete [] pVals;
} /* end Initialize() */

/******************************************************************************
CreateComm()

Creates a model-chromosome communicator.
******************************************************************************/
void ChromosomePool::CreateComm(ModelABC * pModel)
{
   NEW_PRINT("ModelChromoComm", 1);
   m_pComm  = new ModelChromoComm(pModel);   
   MEM_CHECK(m_pComm);
} /* end CreateComm() */

/******************************************************************************
ConvertChromosome()

Utilizes the model-chromosome communicator to convert a chromosome into it's
equivalent parameter group, allowing for more user-friendly output.
******************************************************************************/
ParameterGroup * ChromosomePool::ConvertChromosome(Chromosome * pChromo)
{
   return(m_pComm->ConvertChromosome(pChromo));
}/* end ConvertChromosome() */

/******************************************************************************
WriteMetrics()

Write out setup and metrics for the pool.
******************************************************************************/
void ChromosomePool::WriteMetrics(FILE * pFile)
{
   int i;
   ParameterGroup * pGroup;

   pGroup = ConvertChromosome(GetBestFit());

   fprintf(pFile, "Population Size         : %d\n", m_PoolSize);
   fprintf(pFile, "Number of Elites        : %d\n", m_NumSurvivors);
   fprintf(pFile, "Initialization Method   : ");
   if      (m_InitType == RANDOM_INIT)   { fprintf(pFile, "Random\n");}
   else if (m_InitType == QUAD_TREE_INIT){ fprintf(pFile, "Quad-Tree\n");}
   else if (m_InitType == LHS_INIT)      { fprintf(pFile, "Latin Hypercube Sampling\n");}
   else                                  { fprintf(pFile, "Unknown\n");}

   for(i = 0; i < m_Proto->GetNumGenes(); i++)
   {
      fprintf(pFile, "%-12s Mutations : %d\n", 
      pGroup->GetParamPtr(i)->GetName(), m_pMutCount[i]);
   }/* end for() */   
}/* end WriteMetrics() */
