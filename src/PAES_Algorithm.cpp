/******************************************************************************
File       : PAES_Algorithm.cpp
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)

PAES - Pareto archived evolution strategy

Version History
12-29-17    lsm   created file
******************************************************************************/
#include <mpi.h>
#include <string.h>
#include "PAES_Helpers.h"
#include "PAES_ObjFunc.h"
#include "PAES_Genetics.h"
#include "PAES_Algorithm.h"
#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ObjectiveFunction.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
PAES::PAES(ModelABC * pModel)
{
   int myrank;
   char fname[DEF_STR_SZ];

   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pNonDom = NULL;
   m_pDom = NULL;
   m_precision = NULL;
   m_NumNonDom = 0;
   m_NumDom = 0;
   m_pbuf = NULL;
   m_fbuf = NULL;
   m_pbufsize = 0;
   m_fbufsize = 0;
   
   /* make sure output file name is unique in parallel env */
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(fname, "PAES%d.stdout", myrank);
   m_cout = fopen(fname, "w");

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by PAES and it's member variables.
******************************************************************************/
void PAES::Destroy(void)
{
   ArchiveStruct * pCur, *pDel;

   for (pCur = m_pNonDom; pCur != NULL;)
   {
      pDel = pCur;
      pCur = pCur->pNext;
      delete[] pDel->F;
      delete[] pDel->X;
      delete pDel;
   }

   for (pCur = m_pDom; pCur != NULL;)
   {
      pDel = pCur;
      pCur = pCur->pNext;
      delete[] pDel->F;
      delete[] pDel->X;
      delete pDel;
   }

   delete m_problem;
   delete m_archiveOfSolutions;
   delete m_adaptiveGrid;
   delete [] m_precision;

   delete [] m_pbuf;
   delete [] m_fbuf;

   fclose(m_cout);

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PAES.
******************************************************************************/
void PAES::Calibrate(void)
{
   Optimize();
} /* end Calibrate() */

/******************************************************************************
initialize()
******************************************************************************/
void PAES::initialize(PAES_MultiobjectiveProblem * problemToSolve, PAES_MutationOperator mutationOperator)
{
   int nprocs;
   m_problem = problemToSolve;
   m_mutationOperator = mutationOperator;

   // default values
   m_depth = 4; 
   m_maximumArchiveLength = 150;
   m_numberOfIterations = 1000;
   m_mutationProbability = 0.8; 
   m_seed = 1; 
   m_printFrequency = 100;
   m_distributionIndexForMutation = 20.0; // for polynomial mutation
   m_perturbationForMutation = 10.0; // for uniform mutation
    
   InitFromFile(GetInFileName());
   
   m_numberOfFitnessEvaluations = 0;  

   m_random.initrandom(m_seed);
   m_random.m_Rseed = (float)(m_random.randreal2());
   m_random.randomize();

   m_archiveOfSolutions = new PAES_Population(0, m_maximumArchiveLength, &m_random, m_problem, m_cout);
                                        
   m_adaptiveGrid = new PAES_AdaptiveGrid(m_depth, m_problem->m_numberOfFunctions);

   if (!m_adaptiveGrid)
   {
      fprintf(m_cout, "Paes::Paes-> Not enough memory for the adaptiveGrid object\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */

   //allocate space for parallel communication buffers, if necessary
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   if(nprocs > 1)
   {
      m_pbufsize = m_printFrequency * (m_problem->m_numberOfVariables);
      m_fbufsize = m_printFrequency * (m_problem->m_numberOfFunctions);
      m_pbuf = new double[m_pbufsize];
      m_fbuf = new double[m_fbufsize];
   }
}  /* end PAES::initialize() */

/******************************************************************************
PAES::start()
******************************************************************************/
void PAES::start(void)  
{
   ParameterGroup * pGroup;
   StatusStruct pStatus;
   int innerEvalCount;
   double * X, * F;
   int num, nobj;
   int  result;
   int  mutations;
   long iterations;

   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();
   nobj = m_pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);

   fprintf(m_cout, "START\n" );

   m_startTime = time(NULL);

   /* initial solution */  
   X = new double[num];
   F = new double[nobj];
   WriteInnerEval(WRITE_GA, 1, '.');
   m_currentSolution = new PAES_Individual(m_problem, &m_random, m_cout);
   m_problem->evaluate(m_currentSolution, X, num, F, nobj);
   result = UpdateArchive(X, num, F, nobj);
   if (result == ARCHIVE_NON_DOM)
   {
      WriteInnerEval(1, 1, '+');
   }/* end if() */
   else
   {
      WriteInnerEval(1, 1, '-');
   }/* end else() */
   m_numberOfFitnessEvaluations++;
   addToArchive(new PAES_Individual(m_currentSolution));
   WriteInnerEval(WRITE_ENDED, 0, '.');

   pStatus.pct = 0.00;
   pStatus.numRuns = 1;
   WriteMultiObjRecord(m_pModel, 0, m_pNonDom, pStatus.pct);
   WriteStatus(&pStatus);

   //main optimization loop   
   iterations = 0;
   innerEvalCount = 0;
   while (iterations < m_numberOfIterations)
   {
      if (IsQuit() == true){ break; }

      if (innerEvalCount == 0)
      {
         WriteInnerEval(WRITE_GA, m_printFrequency, '.');
      }

      iterations++;
      pStatus.curIter = m_CurIter = iterations;

      if ((innerEvalCount + 1) == m_printFrequency)
      {
         fprintf(m_cout, "ITERATION = %d\n", (int)iterations);
         fprintf(m_cout, "   Archive Size = : %d \n", m_archiveOfSolutions->getPopulationSize());
         m_currentSolution->printFitness(m_cout);
         fprintf(m_cout, "\n");
      } /* end if() */
      // STEP 1. Copy the current solution
      m_mutantSolution = new PAES_Individual(m_currentSolution);
      if (m_mutantSolution == NULL) 
      {
         fprintf(m_cout, "Paes::start->Error allocating memory for m_mutantSolution\n");
         fclose(m_cout);
         exit(-2);
      } // if
      // STEP 2. Apply at least one mutation
      mutations = 0;
      while(mutations < 1)
      {
         switch (m_mutationOperator)
         {
            case PAES_BIT_FLIP:
               mutations = m_mutantSolution->bitFlipMutation(m_mutationProbability);
               break;
            case PAES_RANDOM:
               mutations = m_mutantSolution->randomMutation(m_mutationProbability);
               break;
            case PAES_POLYNOMIAL:  
               mutations = m_mutantSolution->polynomialMutation(m_mutationProbability, m_distributionIndexForMutation);
               break;
            case PAES_UNIFORM:
               mutations = m_mutantSolution->uniformMutation(m_mutationProbability, m_perturbationForMutation);
               break;
            default:
               fprintf(m_cout, "Paes::start-> mutant operator (%d) undefined\n", m_mutationOperator);
               fclose(m_cout);
               exit(-1);
         } /* end switch */
      }/* end while() */

      if (mutations > 0) 
      {
         X = new double[num];
         F = new double[nobj];

         m_problem->evaluate(m_mutantSolution, X, num, F, nobj);

         result = UpdateArchive(X, num, F, nobj);
         if (result == ARCHIVE_NON_DOM)
         {
            WriteInnerEval(innerEvalCount+1, m_printFrequency, '+');
         }/* end if() */
         else
         {
            WriteInnerEval(innerEvalCount+1, m_printFrequency, '-');
         }/* end else() */

         /* ---------------------------------------------------------
         Test for dominance and update archive accordingly.
         --------------------------------------------------------- */
         testAndUpdate();
      } /* end if() */
      else
      {
         delete m_mutantSolution;
      }
      pStatus.pct = ((float)100.00*(float)(iterations + 1)) / (float)m_numberOfIterations;
      pStatus.numRuns = (iterations + 1);
      WriteStatus(&pStatus);

      innerEvalCount++;
      if(innerEvalCount >= m_printFrequency)
      {
         innerEvalCount = 0;
         WriteInnerEval(WRITE_ENDED, 0, '.');
         WriteMultiObjRecord(m_pModel, (iterations), m_pNonDom, pStatus.pct);
      } /* end if() */      
   } /* end for() */  
   m_endTime = time(NULL);
}  /* end PAES::start() */

/******************************************************************************
PAES::start_parallel()
******************************************************************************/
void PAES::start_parallel(int myrank, int nprocs)
{
   ParameterGroup * pGroup;
   StatusStruct pStatus;
   double * X, * F;
   int num, nobj;
   int  i, j, result;
   int  mutations;
   long iterations;
   double * pMyPrm;
   double * pMyObj;
   int nv;
   int nf;

   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();
   nobj = m_pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);

   fprintf(m_cout, "START\n");

   m_startTime = time(NULL);

   /* initial solution */  
   X = new double[num];
   F = new double[nobj];
   
   if(myrank == 0)
   {
      WriteInnerEval(WRITE_GA, nprocs, '.');
   }
   m_currentSolution = new PAES_Individual(m_problem, &m_random, m_cout);
   m_problem->evaluate(m_currentSolution, X, num, F, nobj);
   delete m_currentSolution;

   /* ------------------------------------------------------------------------
   Each processor has now completed an initial evaluation. Collect the results
   and process them one by one.
   ------------------------------------------------------------------------*/
   nv = m_problem->m_numberOfVariables;
   nf = m_problem->m_numberOfFunctions;

   //each processor packs a different portion of m_pbuf and m_pfbuf
   pMyPrm = &m_pbuf[myrank * nv]; 
   pMyObj = &m_fbuf[myrank * nf];

   for(i = 0; i < nv; i++)
   {
      pMyPrm[i] = X[i];
   } 
   for(i = 0; i < nf; i++)
   {
      pMyObj[i] = F[i];
   } 
   delete [] X;
   delete [] F;

   /* gather the results */
   for(i = 0; i < nprocs; i++)
   {
      pMyPrm = &m_pbuf[i * nv];
      pMyObj = &m_fbuf[i * nf];
      MPI_Bcast(pMyPrm, nv, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(pMyObj, nf, MPI_DOUBLE, i, MPI_COMM_WORLD);
   }/* end for() */
   
   /* process the results */
   for(i = 0; i < nprocs; i++)
   {
      X = new double[num];
      F = new double[nobj];

      for(j = 0; j < nv; j++)
      {
         X[j] = m_pbuf[i*nv + j];
      } 
      for(j = 0; j < nf; j++)
      {
         F[j] = m_fbuf[i*nf + j];
      } 

      result = UpdateArchive(X, num, F, nobj);
      if (result == ARCHIVE_NON_DOM)
      {
         if(myrank == 0)
         {
            WriteInnerEval(1, 1, '+');
         }
      }/* end if() */
      else
      {
         if(myrank == 0)
         {
            WriteInnerEval(1, 1, '-');
         }
      }/* end else() */
      m_numberOfFitnessEvaluations++;
   
      if(i == 0)
      {
         m_currentSolution = new PAES_Individual(X, F, m_problem, &m_random, m_cout);
         addToArchive(new PAES_Individual(m_currentSolution));
      }
      else
      {
         m_mutantSolution = new PAES_Individual(X, F, m_problem, &m_random, m_cout);
         testAndUpdate();
      }
   }/* end for (each processor) */

   if(myrank == 0)
   { 
      WriteInnerEval(WRITE_ENDED, 0, '.');

      pStatus.pct = 0.00;
      pStatus.numRuns = nprocs;
      WriteMultiObjRecord(m_pModel, 0, m_pNonDom, pStatus.pct);
      WriteStatus(&pStatus);
   }/* end if() */

   //main optimization loop
   for(iterations = 0; iterations < m_numberOfIterations; iterations += m_printFrequency)
   {
      if (IsQuit() == true){ break; }

      if (myrank == 0)
      {
         WriteInnerEval(WRITE_GA, m_printFrequency, '.');
      }

      /* parallel runs of the model */
      X = new double[num];
      F = new double[nobj];
      for(j = myrank; j < m_printFrequency; j+=nprocs)
      {
         // STEP 1. Copy the current solution
         m_mutantSolution = new PAES_Individual(m_currentSolution);
         if (m_mutantSolution == NULL)
         {
            fprintf(m_cout, "Paes::start->Error allocating memory for m_mutantSolution\n");
            fclose(m_cout);
            exit(-2);
         } /* end if() */

         // STEP 2. Apply at least one mutation
         mutations = 0;
         while (mutations < 1)
         {
            switch (m_mutationOperator)
            {
               case PAES_BIT_FLIP:
                  mutations = m_mutantSolution->bitFlipMutation(m_mutationProbability);
                  break;
               case PAES_RANDOM:
                  mutations = m_mutantSolution->randomMutation(m_mutationProbability);
                  break;
               case PAES_POLYNOMIAL:
                  mutations = m_mutantSolution->polynomialMutation(m_mutationProbability, m_distributionIndexForMutation);
                  break;
               case PAES_UNIFORM:
                  mutations = m_mutantSolution->uniformMutation(m_mutationProbability, m_perturbationForMutation);
                  break;
               default:
                  fprintf(m_cout, "Paes::start-> mutant operator (%d) undefined\n", m_mutationOperator);
                  fclose(m_cout);
                  exit(-1);
            } /* end switch */
         }/* end while() */

         // Run the model
         m_problem->evaluate(m_mutantSolution, X, num, F, nobj);
         delete m_mutantSolution;

         // save results into shared parallel buffers (m_pbuf and m_fbuf)
         pMyPrm = &m_pbuf[j * nv];
         pMyObj = &m_fbuf[j * nf];
         for (i = 0; i < nv; i++)
         {
            pMyPrm[i] = X[i];
         }
         for (i = 0; i < nf; i++)
         {
            pMyObj[i] = F[i];
         }
      }/* end for() */
      delete[] X;
      delete[] F;

      /* gather the results */
      for (j = 0; j < m_printFrequency; j++)
      {
         pMyPrm = &m_pbuf[j * nv];
         pMyObj = &m_fbuf[j * nf];
         MPI_Bcast(pMyPrm, nv, MPI_DOUBLE, (j%nprocs), MPI_COMM_WORLD);
         MPI_Bcast(pMyObj, nf, MPI_DOUBLE, (j%nprocs), MPI_COMM_WORLD);
      }/* end for() */

      /* process results */
      for (j = 0; j < m_printFrequency; j++)
      {
         X = new double[num];
         F = new double[nobj];

         for (i = 0; i < nv; i++)
         {
            X[i] = m_pbuf[j*nv + i];
         }
         for (i = 0; i < nf; i++)
         {
            F[i] = m_fbuf[j*nf + i];
         }

         result = UpdateArchive(X, num, F, nobj);
         if(myrank == 0)
         {
            if (result == ARCHIVE_NON_DOM)
            {
               WriteInnerEval(j + 1, m_printFrequency, '+');
            }/* end if() */
            else
            {
               WriteInnerEval(j + 1, m_printFrequency, '-');
            }/* end else() */
         }/* end if() */

         /* ---------------------------------------------------------
         Test for dominance and update archive accordingly.
         --------------------------------------------------------- */
         m_mutantSolution = new PAES_Individual(X, F, m_problem, &m_random, m_cout);
         testAndUpdate();
      }/* end for() */      

      pStatus.curIter = m_CurIter = iterations+m_printFrequency;

      fprintf(m_cout, "ITERATION = %d\n", (int)(iterations+m_printFrequency));
      fprintf(m_cout, "   Archive Size = : %d \n", m_archiveOfSolutions->getPopulationSize());
      m_currentSolution->printFitness(m_cout);
      fprintf(m_cout, "\n");

      if(myrank == 0)
      {
         pStatus.pct = ((float)100.00*(float)(iterations + +m_printFrequency)) / (float)m_numberOfIterations;
         pStatus.numRuns = (iterations + nprocs + m_printFrequency);
         WriteStatus(&pStatus);

         WriteInnerEval(WRITE_ENDED, 0, '.');
         WriteMultiObjRecord(m_pModel, (iterations+m_printFrequency), m_pNonDom, pStatus.pct);
      }
   } /* end for() */  

   m_endTime = time(NULL);
}/* end start_parallel() */


/******************************************************************************
PAES::testEncodeDecode()

Perform dominance test and update archive.
******************************************************************************/
void PAES::testEncodeDecode(void)
{
   PAES_BinaryRealGene * theGene = new PAES_BinaryRealGene(m_problem->m_bitsPerVariable[0], m_problem->m_lowerLimit[0], m_problem->m_upperLimit[0], &m_random, m_cout);
   double decoded_true; 
   double decoded_check;
   decoded_true = theGene->getRealAllele();
   theGene->encodeFromReal(decoded_true);
   decoded_check = theGene->getRealAllele();

   fprintf(m_cout, "true: %E\n", decoded_true);
   fprintf(m_cout, "test: %E\n", decoded_check);

}/* end testEncodeDecode() */

/******************************************************************************
PAES::testAndUpdate()

Perform dominance test and update archive.
******************************************************************************/
void PAES::testAndUpdate(void)
{
   int result;
   bool mutantArchived;

   mutantArchived = false;
   m_numberOfFitnessEvaluations++;

   result = m_currentSolution->dominanceTest(m_mutantSolution);
   if (result == -1)
   {
      *m_currentSolution = *m_mutantSolution;
      m_adaptiveGrid->updateGridLocations(m_archiveOfSolutions, m_mutantSolution);
      archiveSolution(m_mutantSolution);
   } /* end if () */
   else if (result == 0) 
   {
      result = compareToArchive(m_mutantSolution);
      // Mutant is not dominated by the archive
      if (result != -1) 
      { 
         m_adaptiveGrid->updateGridLocations(m_archiveOfSolutions, m_mutantSolution);
         mutantArchived = archiveSolution(m_mutantSolution);
         if (mutantArchived && (result == 0)) 
         {
            int location1 = m_adaptiveGrid->findLocation(m_currentSolution);
            int location2 = m_adaptiveGrid->findLocation(m_mutantSolution);
            if (m_adaptiveGrid->m_hypercube[location2] <= m_adaptiveGrid->m_hypercube[location1])
            {
               *m_currentSolution = *m_mutantSolution;
            }
         } /* end if() */
      } /* end if() */
      else
      {
         delete m_mutantSolution;
      }
   } /* end else if() */
   // (result == 1)
   else 
   { 
      delete m_mutantSolution;
   } /* end else() */
}/* end PAES::testAndUpdate() */

/******************************************************************************
PAES::addToArchive()
******************************************************************************/
void PAES::addToArchive(PAES_Individual *solution)
{
   m_archiveOfSolutions->addIndividual(solution);
} /* end PAES::addToArchive() */

/******************************************************************************
PAES::compareToArchive()

Checks if the solution is dominated by any member of the archive.

Returns 0 if the solution is non-dominated, -1 if it is dominated, +1 if it
dominates all the members of the archive.
 ******************************************************************************/
int PAES::compareToArchive(PAES_Individual *solution) 
{
   int counter;
   int result;
  
   counter = 0;
   result  = 0;
  
   while ((counter < m_archiveOfSolutions->getPopulationSize()) && (result != -1)) 
   {
      result = solution->dominanceTest(m_archiveOfSolutions->getIth(counter));
    
      counter ++;
   } /* end while() */
  
   return result;
} /* end PAES::compareToArchive() */

/******************************************************************************
PAES::archiveSolution()
******************************************************************************/
bool PAES::archiveSolution(PAES_Individual *solution) 
{
  int result;
  bool finish;
  int  counter;
  bool storeNewSolution;
  
  finish = false;
  storeNewSolution = true;
  counter = 0;

   /* -------------------------------------------------------------
   CASE 1. The archive is empty
   Action: add solution to the archive
   ------------------------------------------------------------- */
   if (m_archiveOfSolutions->getPopulationSize() == 0)
   {
      addToArchive(solution);
      finish = true;
   } /* end if() */
  
   while ((counter < m_archiveOfSolutions->getPopulationSize()) && !finish) 
   {
      if (solution->identicalFitness(m_archiveOfSolutions->getIth(counter))) 
      {
         // There is a solution with the same fitness vector
         finish = true;
         storeNewSolution = false;
      } // if
      else 
      {
         result = solution->dominanceTest(m_archiveOfSolutions->getIth(counter));
         if (result == -1) 
         {
            finish = true; 
            storeNewSolution = false; // The solution is not added to the archive            
         } /* end if() */
         else if (result == 1) 
         {
            m_archiveOfSolutions->deleteIth(counter);
         } /* end else if() */
         else  
         {
            /* no-op */;
         } /* end else() */
      } /* end else() */
      counter ++;
   } /* end while() */

   if(storeNewSolution) 
   {
      // If the archive is not full, add the solution
      if (m_archiveOfSolutions->getPopulationSize() < m_maximumArchiveLength)
      {
         m_archiveOfSolutions->addIndividual(solution);
      }
      else 
      {
         // The archive is full
         int i;
         int location;
         // If the solution is not in the most crowded region, add it
         location = solution->m_gridLocation;
         if (location == m_adaptiveGrid->m_mostCrowdedHypercube) 
         {
            // archiveOfSolutions_->addIndividual(solution);
            delete solution;
            storeNewSolution = false;
         } // if
         else 
         { // Find and replace an individual of the most crowded region 
            bool finish;
            finish = false;
            i = 0;
            while (!finish) 
            {
               location = (m_archiveOfSolutions->getIth(i))->m_gridLocation;
               if (location == m_adaptiveGrid->m_mostCrowdedHypercube) 
               {
                  m_archiveOfSolutions->deleteIth(i);
                  m_archiveOfSolutions->addIndividual(solution);
                  finish = true;
               } /* end if() */
               i++;
               if (i == m_archiveOfSolutions->getPopulationSize())
               {
                  finish = true;
               }
            } /* end while() */
         } /* end else() */
      } /* end else() */
   } /* end if() */
   else 
   {
      delete solution;
   } /* end else() */
  
   return storeNewSolution;
} /* end PAES::archiveSolution() */

/******************************************************************************
PAES::printToFiles()

Prints the values of the variables and the objective functions to the files 
passed as parameters.

genotypeFileName: The name of the file to store the decision variables
fitnessFileName: The name of the file to store the values of the obj. functions
******************************************************************************/
void PAES::printToFiles(char * genotypeFileName, char * fitnessFileName) 
{
   m_archiveOfSolutions->printGenotype(genotypeFileName);
   m_archiveOfSolutions->printFitness(fitnessFileName);
} /* end PAES::printToFiles() */

/******************************************************************************
UpdateArchive()

Update the dominated and non-dominated archives with latest sample.
******************************************************************************/
int PAES::UpdateArchive(double * pX, int nX, double * pF, int nF)
{
   int i;
   double Fcur, Ftst;
   ArchiveStruct * pArch, *pCur, *pPrev, *pNxt;
   bool bDominates, bIsDominated, bMarkForInsertion;
   pArch = new ArchiveStruct;
   pArch->F = pF;
   pArch->X = pX;
   pArch->nX = nX;
   pArch->nF = nF;
   pArch->pNext = NULL;

   //first entry is always non-dominated
   if ((m_NumDom == 0) && (m_NumNonDom == 0))
   {
      m_pDom = NULL;
      m_pNonDom = pArch;
      m_NumNonDom++;
      return ARCHIVE_NON_DOM;
   }

   //assume solution is non-dominated until we discover otherwise
   bMarkForInsertion = true;

   //compare against current list of non-dominated solutions
   ArchiveStruct * pDummy;
   for (pCur = m_pNonDom; pCur != NULL;)
   {
      //save next item since pCur->pNext may be changed during processing
      pDummy = pCur->pNext;

      //does new solution (Ftst) dominate the existing solution (Fcur)?
      bDominates = true;
      for (i = 0; i < pArch->nF; i++)
      {
         Fcur = pCur->F[i];
         Ftst = pArch->F[i];
         if (Fcur < Ftst)
         {
            bDominates = false;
            break;
         }/* end if() */
      }/* end for() */

      //is new solution (Ftst) dominated by an existing solution (Fcur)?
      if (bDominates == false)
      {
         bIsDominated = true;
         for (i = 0; i < pArch->nF; i++)
         {
            Fcur = pCur->F[i];
            Ftst = pArch->F[i];
            if (Ftst < Fcur)
            {
               bIsDominated = false;
               break;
            }/* end if() */
         }/* end for() */
      }/* end if() */

      /* -----------------------------------------------------------------------------
      Existing solution is dominated. Remove it from list of non-dominated solutions
      and mark new solution for insertion into the non-dominated list.
      ----------------------------------------------------------------------------- */
      if (bDominates == true)
      {
         //solution to be removed is at head of list.
         if (pCur == m_pNonDom)
         {
            m_pNonDom = pCur->pNext;
            m_NumNonDom--;
         }/* end if() */
         else //somewhere in middle or at end.
         {
            for (pPrev = m_pNonDom; pPrev->pNext != pCur;)
            {
               pPrev = pPrev->pNext;
            }
            pNxt = pCur->pNext;
            pPrev->pNext = pNxt; //removes item from list
            m_NumNonDom--;
         }/* end else() */

         //insert at head of dominated list
         pCur->pNext = NULL;
         pNxt = m_pDom;
         m_pDom = pCur;
         m_pDom->pNext = pNxt;
         m_NumDom++;
      }/* end if() */
      /* -----------------------------------------------------------------------------
      New solution is dominated. Make note so that it is not inserted into the list
      of non-dominated solutions.
      ----------------------------------------------------------------------------- */
      else if (bIsDominated == true)
      {
         bMarkForInsertion = false;
      }/* end if() */

      //advance to next item
      pCur = pDummy;
   }/* end for() */

   //insert new solution into list of non-dominated solutions?
   if (bMarkForInsertion == true)
   {
      //insert at head of non-dominated list
      pNxt = m_pNonDom;
      m_pNonDom = pArch;
      m_pNonDom->pNext = pNxt;
      m_NumNonDom++;
      return ARCHIVE_NON_DOM;
   }/* end if() */
   else
   {
      //insert at head of dominated list
      pNxt = m_pDom;
      m_pDom = pArch;
      m_pDom->pNext = pNxt;
      m_NumDom++;
      return ARCHIVE_DOM;
   }/* end else() */
}/* end PAES::UpdateArchive() */

/******************************************************************************
PAES::InitFromFile()

Reads configurations parameters from a file
******************************************************************************/
void PAES::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   int num, nprocs;   

   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   num = m_pModel->GetParamGroupPtr()->GetNumParams();

   m_depth = 4;
   m_maximumArchiveLength = 150;
   m_numberOfIterations = 1000;
   m_seed = GetRandomSeed();

   m_mutationProbability = 0.8;

   m_printFrequency = 100*nprocs;
   m_distributionIndexForMutation = 20.0; // for polynomial mutation
   m_perturbationForMutation = 10.0; // for uniform mutation

   //read in PAES configuration
   pFile = fopen(pFileName, "r");
   if (pFile == NULL)
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open PAES config. file. Using Defaults");
      return;
   }/* end if() */ 

   //make sure correct tokens are present
   if (CheckToken(pFile, "BeginPAES", pFileName) == true)
   {
      FindToken(pFile, "EndPAES", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginPAES", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while (strstr(line, "EndPAES") == NULL)
      {
         if (strstr(line, "Depth") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_depth);
         }/*end else if() */
         else if (strstr(line, "MaxArchiveLength") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_maximumArchiveLength);
         }/*end else if() */
         else if (strstr(line, "NumberOfIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_numberOfIterations);
         }/*end else if() */
         else if (strstr(line, "MutationProb") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_mutationProbability);
         }/*end else if() */
         else if (strstr(line, "PrintFreq") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_printFrequency);
         }/*end else if() */
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */

   /* avoid wasting processors */
   if(m_printFrequency < nprocs)
   {
      fprintf(m_cout, "PAES::InitFromFile() -> Increasing print frequency from %d to %d to avoid idling processors\n", m_printFrequency, nprocs);
      m_printFrequency = nprocs;
   }

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
PAES::printStatistics()

Print execution statistics
******************************************************************************/
void PAES::printStatistics(void) 
{
   fprintf(m_cout, "   RESULTS\n");
   fprintf(m_cout, "-------------\n");
   fprintf(m_cout, "Time: %f seconds\n", (double)((double)m_endTime - (double)m_startTime));
   fprintf(m_cout, "Evaluations: %d\n", m_numberOfFitnessEvaluations);
   fprintf(m_cout, "\nEND\n");;
} /* end PAES::printStatistics() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void PAES::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Pareto Archived Evolution Strategy\n");
   fprintf(pFile, "Depth                   : %d\n", m_depth);
   fprintf(pFile, "Max Archive Length      : %d\n", m_maximumArchiveLength);
   fprintf(pFile, "Number of Iterations    : %d\n", m_numberOfIterations);
   fprintf(pFile, "Mutation Probability    : %f\n", m_mutationProbability);
   fprintf(pFile, "Print Frequency         : %d\n", m_printFrequency);
   fprintf(pFile, "Num Fitness Evals       : %d\n", m_numberOfFitnessEvaluations);
   fprintf(pFile, "Non-Dominated Solutions : %d\n", m_NumNonDom);
   fprintf(pFile, "Dominated Solutions     : %d\n", m_NumDom);
  
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
PAES::Optimize()
******************************************************************************/
void PAES::Optimize(void)
{
   int myrank, nprocs;
   StatusStruct pStatus;
   PAES_MultiobjectiveProblem * problemToSolve;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   problemToSolve = new PAES_OstrichIf(PAES_BINARY_REAL, m_pModel, m_cout);

   initialize(problemToSolve, PAES_BIT_FLIP);

   if(myrank == 0)
   {
      WriteSetup(m_pModel, "PAES - Pareto Archived Evolution Strategy");
      //write banner
      WriteBanner(m_pModel, "gen   ", "Precent Complete");
   }

   if(nprocs < 2)
   {
      start();
   }
   else
   {
      start_parallel(myrank, nprocs);
   }

   if(myrank == 0)
   {
      printStatistics();

      printToFiles((char *)"PAES_Parameters.out", (char *)"PAES_ObjFunctions.out");
 
      WriteMultiObjOptimal(m_pModel, m_pNonDom, m_pDom);

      /* final status */
      pStatus.curIter = m_CurIter;
      pStatus.maxIter = m_numberOfIterations;
      pStatus.pct = ((float)100.00*(float)(m_CurIter + 1)) / (float)m_numberOfIterations;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //write algorithm metrics
      WriteAlgMetrics(this);
   }/* end if() */
} /* end Optimize() */

/******************************************************************************
PAES_Program()

Calibrate or optimize the model using PAES.
******************************************************************************/
void PAES_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("PAES", 1);
   PAES * TheAlg = new PAES(model);
   MEM_CHECK(TheAlg);

   if (model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end PAES_Program() */

