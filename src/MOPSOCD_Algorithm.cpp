/******************************************************************************
File     : MOPSOCD.h

Authors  : L. Shawn Matott, C.R. Raquel and P.C. Naval, Jr.

Copyrights: 2005 C.R. Raquel and P.C. Naval, Jr., University of the Philippines-Diliman
--- original implementation as standalone C program
--- MOPSO-CD Version 0.5b Feb. 17, 2006
2017, L. Shawn Matott
--- adaptation into OSTRICH toolkit

Multi-Objective Particle Swarm Optimization with Crowding Distance (MOPSOCD).

This is an implementation of MOPSO-CD,a multiobjective particle swarm
optimization algorithm using crowding distance.

For details please see:
C.R. Raquel and P.C. Naval, Jr., An Effective Use of Crowding
Distance in Multiobjective Particle Swarm Optimization.
In Proc. of Genetic and Evolutionary Computation Conference
(GECCO 2005), Washington DC, June 2005.

E-mail address       : cvmig@engg.upd.edu.ph
Version              : mopsocd05b
Last updated         : Fri Feb 17, 2006

Random Generator Source code has been taken from Random Library found
at http://www.swin.edu.au/astronomy/pbourke/software/random/

Permission to use MOPSO-CD codes is hereby granted for academic and
research purposes ONLY. Commercial usage of these codes is prohibited
without prior knowledge of the authors.  In no way will the authors
be held responsible for any possible faulty operation of
software/hardware arising from the use of these codes.

Version History
12-28-17    lsm   added copyright information and initial comments.
******************************************************************************/

#include <mpi.h>
#include "MOPSOCD_Algorithm.h"
#include "Model.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
MOPSOCD::MOPSOCD(ModelABC * pModel)
{
   ParameterGroup * pGroup;

   pGroup = pModel->GetParamGroupPtr();
   m_NumVar = pGroup->GetNumParams();
   m_NumFun = pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);
   
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pNonDomList = NULL;
   m_pDomList = NULL;
   m_NonDomListSize = 0;
   m_DomListSize = 0;

   m_NumNonDom = 0;
   m_PopSize = 0;
   m_MaxGen = 0;
   m_CurIter = 0;
   m_ArchiveSize = 0;
   
   m_PI = 4.0*atan(1.0);
   m_ProbMut = 0.5;

   m_ArchiveVar = NULL;
   m_ArchiveFit = NULL;
   m_PopVar = NULL;
   m_PopFit = NULL;
   m_PbestsVar = NULL;
   m_PbestsFit = NULL;
   m_Velocity = NULL;
   m_CrowdDist = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by MOPSOCD and it's member variables.
******************************************************************************/
void MOPSOCD::Destroy(void)
{
   int i;

   ArchiveStruct * pCur, *pDel;

   for (pCur = m_pNonDomList; pCur != NULL;)
   {
      pDel = pCur;
      pCur = pCur->pNext;
      delete[] pDel->F;
      delete[] pDel->X;
      delete pDel;
   }

   for (pCur = m_pDomList; pCur != NULL;)
   {
      pDel = pCur;
      pCur = pCur->pNext;
      delete[] pDel->F;
      delete[] pDel->X;
      delete pDel;
   }

   for (i = 0; i < m_ArchiveSize; i++)
   {
      delete [] m_ArchiveVar[i];
      delete [] m_ArchiveFit[i];
   }
   delete[] m_ArchiveVar;
   delete[] m_ArchiveFit;

   for (i = 0; i < m_PopSize; i++)
   {
      delete[] m_PopVar[i];
      delete[] m_PopFit[i];
   }
   delete[] m_PopVar;
   delete[] m_PopFit;

   for (i = 0; i < m_PopSize; i++)
   {
      delete[] m_PbestsVar[i];
      delete[] m_PbestsFit[i];
   }
   delete[] m_PbestsVar;
   delete[] m_PbestsFit;

   for (i = 0; i < m_PopSize; i++)
   {
      delete[] m_Velocity[i];
   }
   delete[] m_Velocity;

   delete [] m_CrowdDist;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using MOPSOCD.
******************************************************************************/
void MOPSOCD::Calibrate(void)
{
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Search for the pareto front using MOPSOCD.
******************************************************************************/
void MOPSOCD::Optimize(void)
{
   char archiveName[20];
   int i, j, g, myrank, nprocs;
   clock_t  startTime, endTime;
   double duration, clocktime;
   FILE *outfile, *plotfile;

   StatusStruct pStatus;
   ParameterGroup * pGroup;

   InitFromFile(GetInFileName());

   pGroup = m_pModel->GetParamGroupPtr();

   /* get rank and num procs */
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   sprintf(archiveName, "MOPSOCD_archive.out");

   if(myrank == 0)
   {
      WriteSetup(m_pModel, "MOPSOCD - Multi-Objective Particle Swarm Optimization with Crowding Distance");
      //write banner
      WriteBanner(m_pModel, "gen   ", "Percent Complete");   
   
      outfile = fopen("MOPSOCD_output.out", "w");
      plotfile = fopen("MOPSOCD_plot.out", "w");
   }

   /* Initialize random number generator */
   initialize_rand();

   startTime = clock();

   /* Initialize generation counter */
   m_CurIter = g = 0;

   /* Initialize population with random values */
   initialize_pop();

   /* Initialize velocity */
   initialize_vel();

   //Update run record to indicate onset of PSO sampling 
   if(myrank == 0)
   {
      WriteInnerEval(WRITE_PSO, m_PopSize, '.');
   }

   /* Evaluate particles in population, possibly in parallel */
   if(nprocs > 1)
   {
      evaluate_parallel(myrank, nprocs);
   }
   else
   {
      evaluate();
   }

   if(myrank == 0)
   {
      /* update OSTRICH run record and status file */
      WriteInnerEval(WRITE_ENDED, 0, '.');
      pStatus.pct = ((float)100.00*(float)(g)) / (float)m_MaxGen;
      pStatus.numRuns = (g + 1)*m_PopSize;
      WriteMultiObjRecord(m_pModel, g, m_pNonDomList, pStatus.pct);
      WriteStatus(&pStatus);
   }

   /* Store initial personal bests (both variable and fitness values) of particles */
   store_pbests();

   /* Insert nondominated particles in population into the archive */
   insert_nondom();

   /******* MAIN OPTIMIZATION LOOP *********/
   g = 1;
   while (g <= m_MaxGen)
   {
      clocktime = (clock() - startTime) / (double)CLOCKS_PER_SEC;

      pStatus.curIter = m_CurIter = g;

      if (IsQuit() == true){ break; }

      /* report progress to stdout */
      if ((MOPSOCD_VERBOSE > 0) && (((g % MOPSOCD_PRINT_EVERY) == 0) || (g == m_MaxGen)))
      { 
         if(myrank == 0)
         {
            fprintf(stdout, "Generation: %d      Time: %.2f sec\n", g, clocktime);
            fflush(stdout);
         }
      }

      /* report progress to out file */
      if ((g % MOPSOCD_PRINT_EVERY) == 0 || g == m_MaxGen)
      {
         if(myrank == 0)
         {
            fprintf(outfile, "Generation: %d     Time: %.2f sec\n", g, clocktime);
         }
      }

      /* ------------------------------------------------------------
         Compute crowding distance of each particle in the archive. 
         Must be at least 3 particles in the archive.
      ------------------------------------------------------------ */
      if (m_NumNonDom > 2)
      {
         crowding();
      }

      /* Compute new velocity of each particle in the population */
      compute_velocity();

      /* Maintain particles in the population within the search space */
      maintain_particles();

      /* Mutate particles in the population */
      if (g < (m_MaxGen * m_ProbMut))
      {
         mutate(g);
      }

      //Update run record to indicate onset of PSO sampling 
      if(myrank == 0)
      {
         WriteInnerEval(WRITE_PSO, m_PopSize, '.');
      }

      /* Evaluate particles in the population */
      if(nprocs > 1)
      {
         evaluate_parallel(myrank, nprocs);
      }
      else
      {
         evaluate();
      }

      /* Insert new nondominated particles in pop into archive */
      update_archive();

      /* Update personal bests of particles in the population */
      update_pbests();

      /* write out best so far */
      if (((g % MOPSOCD_PRINT_EVERY) == 0) || (g == m_MaxGen))
      {
         if(myrank == 0)
         {
            fprintf(outfile, "Size of Pareto Set: %d\n\n", m_NumNonDom);
            for (i = 0; i < m_NumNonDom; i++)
            {
               fprintf(outfile, "Parameter Values:\n");
               for (j = 0; j < m_NumVar; j++)
               {
                  fprintf(outfile, "   %E\n", m_ArchiveVar[i][j]);
               }

               fprintf(outfile, "Obj. Function Values:\n");
               for (j = 0; j < m_NumFun; j++)
               {
                  fprintf(outfile, "   %E\n", m_ArchiveFit[i][j]);
               }
               fprintf(outfile, "\n");
            }/* end for() */

            fflush(outfile);
         }/* end if() */
      }/* end if() */

      /* write out GNU plot file */
      if(g == m_MaxGen) 
      {
         if(myrank == 0)
         {
            fprintf(plotfile, "# GNU Plot\n");
            for (i = 0; i < m_NumNonDom; i++) 
            {
               for (j = 0; j < m_NumFun; j++)
               {
                  fprintf(plotfile, "%.4f ", m_ArchiveFit[i][j]);
               }
               fprintf(plotfile, "\n");
            }
            fflush(plotfile);
         }/* end if() */
      }/* end if() */

      /* update OSTRICH run record and status file */
      if(myrank == 0)
      {
         WriteInnerEval(WRITE_ENDED, 0, '.');     
         pStatus.pct = ((float)100.00*(float)(g)) / (float)m_MaxGen;
         pStatus.numRuns = (g + 1)*m_PopSize;
         WriteMultiObjRecord(m_pModel, g, m_pNonDomList, pStatus.pct);
         WriteStatus(&pStatus);
      }

      /* Increment generation counter */
      g++;
   }/* end while() */

   /* Write results to file */
   if(myrank == 0)
   {
      save_results(archiveName);
   }

   /* Compute time duration */
   endTime = clock();
   duration = (endTime - startTime) / (double)CLOCKS_PER_SEC;
   if (myrank == 0)
   {
      fprintf(stdout, "%lf sec\n", duration);
   
      fclose(outfile);
      fclose(plotfile);

      WriteMultiObjOptimal(m_pModel, m_pNonDomList, m_pDomList);
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      //write algorithm metrics
      WriteAlgMetrics(this);
   }
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void MOPSOCD::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : MOPSOCD - Multi-Objective Particle Swarm Optimization with Crowding Distance\n");
   fprintf(pFile, "Max Generations         : %d\n", m_MaxGen);
   fprintf(pFile, "Population Size         : %d\n", m_PopSize);
   fprintf(pFile, "Mutation Rate           : %0.2f\n", m_ProbMut);
   fprintf(pFile, "Archived Solutions      : %d\n", m_NumNonDom);
   fprintf(pFile, "Non-Dominated Solutions : %d\n", m_NonDomListSize);
   fprintf(pFile, "Dominated Solutions     : %d\n", m_DomListSize);

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void MOPSOCD::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   int i;

   m_PopSize = 20;
   m_MaxGen = 50;
   m_ArchiveSize = 1000;
   m_ProbMut = 0.5;

   //read in MOPSOCD configuration
   pFile = fopen(pFileName, "r");
   if (pFile == NULL)
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open MOPSOCD config. file. Using Defaults");
      return;
   }/* end if() */

   //make sure correct tokens are present
   if (CheckToken(pFile, "BeginMOPSOCD", pFileName) == true)
   {
      FindToken(pFile, "EndMOPSOCD", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginMOPSOCD", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while (strstr(line, "EndMOPSOCD") == NULL)
      {
         if (strstr(line, "PopSize") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_PopSize);
         }/*end else if() */
         else if (strstr(line, "MaxGen") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxGen);
         }
         else if (strstr(line, "ArchiveSize") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_ArchiveSize);
         }
         else if (strstr(line, "MutationRate") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_ProbMut);
         }
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */

   // allocate memory
   m_ArchiveVar = new double*[m_ArchiveSize];
   for(i = 0; i < m_ArchiveSize; i++)
   {
      m_ArchiveVar[i] = new double[m_NumVar];
   }

   m_ArchiveFit = new double*[m_ArchiveSize];
   for (i = 0; i < m_ArchiveSize; i++)
   {
      m_ArchiveFit[i] = new double[m_NumFun];
   }

   m_CrowdDist = new double [m_ArchiveSize];

   m_PopVar = new double*[m_PopSize];
   for (i = 0; i < m_PopSize; i++)
   {
      m_PopVar[i] = new double[m_NumVar];
   }

   m_PopFit = new double*[m_PopSize];
   for (i = 0; i < m_PopSize; i++)
   {
      m_PopFit[i] = new double[m_NumFun];
   }

   m_PbestsVar  = new double*[m_PopSize];
   for (i = 0; i < m_PopSize; i++)
   {
      m_PbestsVar[i] = new double[m_NumVar];
   }

   m_PbestsFit  = new double*[m_PopSize];
   for (i = 0; i < m_PopSize; i++)
   {
      m_PbestsFit[i] = new double[m_NumFun];
   }

   m_Velocity = new double*[m_PopSize];
   for (i = 0; i < m_PopSize; i++)
   {
      m_Velocity[i] = new double[m_NumVar];
   }

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
UpdateLists()

Update the dominated and non-dominated lists with latest sample.
******************************************************************************/
int MOPSOCD::UpdateLists(double * pX, int nX, double * pF, int nF)
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
   if ((m_DomListSize == 0) && (m_NonDomListSize == 0))
   {
      m_pDomList = NULL;
      m_pNonDomList = pArch;
      m_NonDomListSize++;
      return ARCHIVE_NON_DOM;
   }

   //assume solution is non-dominated until we discover otherwise
   bMarkForInsertion = true;

   //compare against current list of non-dominated solutions
   ArchiveStruct * pDummy;
   for (pCur = m_pNonDomList; pCur != NULL;)
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
         if (pCur == m_pNonDomList)
         {
            m_pNonDomList = pCur->pNext;
            m_NonDomListSize--;
         }/* end if() */
         else //somewhere in middle or at end.
         {
            for (pPrev = m_pNonDomList; pPrev->pNext != pCur;)
            {
               pPrev = pPrev->pNext;
            }
            pNxt = pCur->pNext;
            pPrev->pNext = pNxt; //removes item from list
            m_NonDomListSize--;
         }/* end else() */

         //insert at head of dominated list
         pCur->pNext = NULL;
         pNxt = m_pDomList;
         m_pDomList = pCur;
         m_pDomList->pNext = pNxt;
         m_DomListSize++;
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
      pNxt = m_pNonDomList;
      m_pNonDomList = pArch;
      m_pNonDomList->pNext = pNxt;
      m_NonDomListSize++;
      return ARCHIVE_NON_DOM;
   }/* end if() */
   else
   {
      //insert at head of dominated list
      pNxt = m_pDomList;
      m_pDomList = pArch;
      m_pDomList->pNext = pNxt;
      m_DomListSize++;
      return ARCHIVE_DOM;
   }/* end else() */
}/* end UpdateLists() */

/******************************************************************************
initialize_rand()

Initialize random number generator from MOPSOCD_Randomlib.h
******************************************************************************/
void MOPSOCD::initialize_rand(void)
{
   unsigned int i, j;

   srand((unsigned int)time((time_t *)NULL));
   i = (unsigned int)(31329.0 * rand() / (RAND_MAX + 1.0));
   j = (unsigned int)(30082.0 * rand() / (RAND_MAX + 1.0));

   MOPSOCD_RandomInitialise(i, j);
}/* end initialize_rand() */

/******************************************************************************
initialize_pop()

Initialize population variables.
******************************************************************************/
void MOPSOCD::initialize_pop(void)
{
   int i, j;
   double xmin, xmax;

   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         xmin = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetLwrBnd();
         xmax = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetUprBnd();
         m_PopVar[i][j] = MOPSOCD_RandomDouble(xmin, xmax);
      }
   }/* end for() */
}/* end initialize_pop() */

/******************************************************************************
initialize_vel()

Initialize population velocity to zero.
******************************************************************************/
void MOPSOCD::initialize_vel(void) 
{
   int i, j;

   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         m_Velocity[i][j] = 0.0;
      }
   }
}/* end initialize_vel() */

/******************************************************************************
evaluate()

Evaluate particles in population.
******************************************************************************/
void MOPSOCD::evaluate(void) 
{
   int result;
   int i, j;
   double * F, * X;
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();

   for (i = 0; i < m_PopSize; i++)
   {
      X = new double[m_NumVar];
      F = new double[m_NumFun];

      for(j = 0; j < m_NumVar; j++)
      {
         X[j] = m_PopVar[i][j];
      }
      pGroup->WriteParams(X);
      m_pModel->Execute(F, m_NumFun);

      for (j = 0; j < m_NumFun; j++)
      {
         m_PopFit[i][j] = F[j];
      }

      result = UpdateLists(X, m_NumVar, F, m_NumFun);
      if (result == ARCHIVE_NON_DOM)
      {
         WriteInnerEval(i + 1, m_PopSize, '+');
      }/* end if() */
      else
      {
         WriteInnerEval(i + 1, m_PopSize, '-');
      }/* end else() */
   }/* end for() */
}/* end evaluate() */

/******************************************************************************
evaluate_parallel()

Evaluate particles in population using parallel processing.
******************************************************************************/
void MOPSOCD::evaluate_parallel(int myrank, int nprocs)
{
   int result;
   int i, j;
   double * F, *X;
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();

   broadcast_population(myrank, nprocs);

   X = new double[m_NumVar];
   F = new double[m_NumFun];
   for (i = myrank; i < m_PopSize; i+=nprocs)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         X[j] = m_PopVar[i][j];
      }
      pGroup->WriteParams(X);
      m_pModel->Execute(F, m_NumFun);

      for (j = 0; j < m_NumFun; j++)
      {
         m_PopFit[i][j] = F[j];
      }
   }/* end for() */
   delete[] X;
   delete[] F;

   gather_results(myrank, nprocs);

   if(myrank == 0)
   {
      for (i = 0; i < m_PopSize; i++)
      {
         X = new double[m_NumVar];
         F = new double[m_NumFun];

         for (j = 0; j < m_NumVar; j++)
         {
            X[j] = m_PopVar[i][j];
         }

         for (j = 0; j < m_NumFun; j++)
         {
            F[j] = m_PopFit[i][j];
         }

         result = UpdateLists(X, m_NumVar, F, m_NumFun);
         if (result == ARCHIVE_NON_DOM)
         {
            WriteInnerEval(i + 1, m_PopSize, '+');
         }/* end if() */
         else
         {
            WriteInnerEval(i + 1, m_PopSize, '-');
         }/* end else() */
      }/* end for() */
   } /* end if() */
}/* end evaluate_parallel() */

/******************************************************************************
broadcast_population()

Use MPI communications to ensure that all processor have the same population.
******************************************************************************/
void MOPSOCD::broadcast_population(int myrank, int nprocs)
{
   double * X;
   int i, j, k, nx;

   nx = m_NumVar*m_PopSize;
   X = new double[nx];

   /* flatten array for easier tranmission */
   k = 0;
   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         X[k] = m_PopVar[i][j];
         k++;
      }
   }/* end for() */

   /* boradcast out supervisor data */
   MPI_Bcast(X, nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   /* unpack the array */
   k = 0;
   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         m_PopVar[i][j] = X[k];
         k++;
      }
   }/* end for() */

   delete [] X;
} /* end broadcast_population() */

/******************************************************************************
gather_results()

Use MPI communications to collect results.
******************************************************************************/
void MOPSOCD::gather_results(int myrank, int nprocs)
{
   int i, j, k, p, kstep;
   int nf;
   double * F; // a flattened array for easier tranmission

   nf = m_NumFun * m_PopSize;
   kstep = m_NumFun * (nprocs - 1);
   F = new double[nf];

   /* share data among processors */
   for(p = 0; p < nprocs; p++)
   {
      /* pack F with processor data */
      if(myrank == p)
      {         
         k = p*m_NumFun;
         for (i = p; i < m_PopSize; i += nprocs)
         {
            for (j = 0; j < m_NumFun; j++)
            {
               F[k] = m_PopFit[i][j];
               k++;
            }/* end for() */
            k += kstep;
         }/* end for() */
      }/* end if() */

      /* broadcast data */
      MPI_Bcast(F, nf, MPI_DOUBLE, p, MPI_COMM_WORLD);
      
      /* upack processor results */
      k = p*m_NumFun;
      for(i = p; i < m_PopSize; i+=nprocs)
      {
         for (j = 0; j < m_NumFun; j++)
         {
            m_PopFit[i][j] = F[k];
            k++;
         }/* end for() */
         k += kstep;
      }/* end for() */
   }/* end for() */

   delete[] F;
} /* end gather_results() */

/******************************************************************************
store_pbests()

Store personal bests (both variable and fitness values) of particles.
******************************************************************************/
void MOPSOCD::store_pbests(void) 
{
   int i, j;

   /* Store variable values of personal bests */
   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         m_PbestsVar[i][j] = m_PopVar[i][j];
      }
   }

   /* Store fitness values of personal bests */
   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumFun; j++)
      {
         m_PbestsFit[i][j] = m_PopFit[i][j];
      }
   }
}/* end store_pbests() */

/******************************************************************************
insert_nondom()

Insert nondominated particles in population into the archive.
******************************************************************************/
void MOPSOCD::insert_nondom()
{
   int i, j, k, total, insertFlag, bottom;
   double * archiveCons, *popCons;

   archiveCons = new double[m_NumVar];
   popCons = new double[m_NumVar];

   for (i = 0; i< m_PopSize; i++)
   {
      insertFlag = 1;

      /* if archive is empty */
      if (m_NumNonDom == 0)
      { 

         /* Insert particle in pop into archive */
         for (j = 0; j < m_NumVar; j++)
         {
            m_ArchiveVar[m_NumNonDom][j] = m_PopVar[i][j];
         }

         for (j = 0; j < m_NumFun; j++)
         {
            m_ArchiveFit[m_NumNonDom][j] = m_PopFit[i][j];
         }

         m_NumNonDom += 1;
      }/* end if() */
      /* if archive is not empty */
      else
      { 
         insertFlag = 1;

         /*for each particle in archive	*/
         for (k = 0; k < m_NumNonDom; k++)
         {
            /* First, check for feasibility */
            for (j = 0; j < m_NumVar; j++)
            {
               popCons[j] = m_PopVar[i][j];
               archiveCons[j] = m_ArchiveVar[k][j];
            }

            /* If both particles are infeasible */
            if ((check_constraints(archiveCons) > 0) && (check_constraints(popCons) > 0))
            {
               delete_particle(k);           /* Delete particle in archive           */
               insertFlag = 0;		/* Do not insert particle in pop        */
               break;
            }
            /* If particle in pop is infeasible */
            else if (check_constraints(popCons) > 0)
            { 
               insertFlag = 0;			   /* Do not insert particle in pop    */
               break;
            }
            /* If particle in archive is infeasible */
            else if (check_constraints(archiveCons) > 0)
            { 
               delete_particle(k); /* Delete particle in archive */

               if ((m_NumNonDom != 0) || (k != (m_NumNonDom - 1)))
               {
                  k--;
               }

               continue;
            }/* end else if() */

            /* Second, check for domination */
            total = 0;

            /* If both are feasible, check for nondomination */
            for (j = 0; j < m_NumFun; j++)
            {
               if (((m_PopFit[i][j] < m_ArchiveFit[k][j]) && (MOPSOCD_OPTIMIZATION_TYPE == 0))
                  || (m_PopFit[i][j] > m_ArchiveFit[k][j]) && (MOPSOCD_OPTIMIZATION_TYPE == 1))
               {
                  total += 1;
               }
            }

            /* If particle in pop dominates     */
            if (total == m_NumFun)   
            {
               delete_particle(k); /* Delete particle in archive       */
            }
            /* If particle in archive dominates */
            else if (total == 0)
            {  
               insertFlag = 0;     /* Do not insert particle in pop    */
               break;
            }
         } /* end for() --- Finished comparing one particle in pop with particles in archive */
      }/* end else() */

      /* Insert particle in pop if it is feasible and nondominated */
      if (insertFlag == 1)
      {
         /* If memory is not yet full, insert particle */
         if (m_NumNonDom < m_ArchiveSize)
         {
            for (j = 0; j < m_NumVar; j++)
            {
               m_ArchiveVar[m_NumNonDom][j] = m_PopVar[i][j];
            }
            for (j = 0; j < m_NumFun; j++)
            {
               m_ArchiveFit[m_NumNonDom][j] = m_PopFit[i][j];
            }
            m_NumNonDom += 1;
         }/* end if() */
         /* If memory is full, select particle to replace */
         else
         {      
            /* Compute crowding distance values of particles in archive*/
            crowding();

            bottom = (unsigned int)((m_NumNonDom - 1) * 0.90);

            /* Randomly select which particle from the most areas to replace */
            k = MOPSOCD_RandomInt(bottom, m_NumNonDom - 1);

            /* Insert new particle into archive */
            for (j = 0; j < m_NumVar; j++)
            {
               m_ArchiveVar[k][j] = m_PopVar[i][j];
            }

            for (j = 0; j < m_NumFun; j++)
            {
               m_ArchiveFit[k][j] = m_PopFit[i][j];
            }
         }/* else() */
      }/* end if() */
   } /* end for() --- Finished comparing particles in pop with particles in archive */

   delete [] archiveCons;
   delete [] popCons;

}/* end insert_nondom() */

/******************************************************************************
delete_particle()

Delete particle in archive.
******************************************************************************/
void MOPSOCD::delete_particle(unsigned int k)
{
   int j;

   /* if infeasible particle is the last in archive or only one particle in archive */
   if ((m_NumNonDom == 1) || (k == (m_NumNonDom - 1))) 
   {
      m_NumNonDom -= 1;
   }
   /* move last particle in archive in place of infeasible particle	*/
   else 
   {       
      for (j = 0; j < m_NumVar; j++)
      {
         m_ArchiveVar[k][j] = m_ArchiveVar[m_NumNonDom - 1][j];
      }

      for (j = 0; j < m_NumFun; j++)
      {
         m_ArchiveFit[k][j] = m_ArchiveFit[m_NumNonDom - 1][j];
      }

      m_NumNonDom -= 1;
   }
}/* end delete_particle() */

/******************************************************************************
check_constraints()

Check for constraint to determine feasibility
******************************************************************************/
unsigned int MOPSOCD::check_constraints(double *consVar)
{
   return 0;
}/* end check_constraints() */

/******************************************************************************
crodwing()

Computes crowding distance values of particles in archive
******************************************************************************/
void MOPSOCD::crowding(void)
{
   int i, f, begin;

   /* initialize crowding distance values */
   for (i = 0; i <= m_NumNonDom; i++)
   {
      m_CrowdDist[i] = 0;
   }

   for (f = 0; f < m_NumFun; f++)
   {
      begin = 0;
      /* Sort fitness values */
      qsortFitness(f, begin, m_NumNonDom);

      /* Compute crowding distance */
      compute_distance(f);
   }

   /* Sort crowding distance values */
   begin = 0;
   qsortCrowd(begin, m_NumNonDom);
}/* end crodwing() */

/******************************************************************************
qsortFitness()

Sort fitness values of particles in archive
******************************************************************************/
void MOPSOCD::qsortFitness(unsigned int f, unsigned int begin, unsigned int lastPart)
{
   unsigned int left = begin + 1;
   unsigned int right = lastPart;
   double pivot = m_ArchiveFit[begin][f];

   int k;
   double *tempF, * tempP, temp;

   tempF = new double[m_NumFun];
   tempP = new double[m_NumVar];

   while (left < right)
   {
      if (m_ArchiveFit[left][f] <= pivot) 
      {
         left++;
      }
      else 
      {
         right--;

         /* swap */
         /* Exchange fitness positions of two particles in the archiveFit */
         for (k = 0; k < m_NumFun; k++)
         {
            tempF[k] = m_ArchiveFit[left][k];
            m_ArchiveFit[left][k] = m_ArchiveFit[right][k];
         }

         for (k = 0; k < m_NumFun; k++)
         {
            m_ArchiveFit[right][k] = tempF[k];
         }

         /* Also exchange particle positions in archiveVar */
         for (k = 0; k < m_NumVar; k++)
         {
            tempP[k] = m_ArchiveVar[left][k];
            m_ArchiveVar[left][k] = m_ArchiveVar[right][k];
         }

         for (k = 0; k < m_NumVar; k++)
         {
            m_ArchiveVar[right][k] = tempP[k];
         }

         /* Also exchange their crowding distance */
         temp = m_CrowdDist[left];
         m_CrowdDist[left] = m_CrowdDist[right];
         m_CrowdDist[right] = temp;
      }
   }

   left--;

   /* Exchange fitness positions of two particles in the array archiveFit */
   for (k = 0; k < m_NumFun; k++)
   {
      tempF[k] = m_ArchiveFit[begin][k];
      m_ArchiveFit[begin][k] = m_ArchiveFit[left][k];
   }

   for (k = 0; k < m_NumFun; k++)
   {
      m_ArchiveFit[left][k] = tempF[k];
   }

   /* Also exchange particle positions in archiveVar */
   for (k = 0; k < m_NumVar; k++)
   {
      tempP[k] = m_ArchiveVar[begin][k];
      m_ArchiveVar[begin][k] = m_ArchiveVar[left][k];
   }

   for (k = 0; k < m_NumVar; k++)
   {
      m_ArchiveVar[left][k] = tempP[k];
   }

   /* Also exchange their crowding distance */
   temp = m_CrowdDist[begin];
   m_CrowdDist[begin] = m_CrowdDist[left];
   m_CrowdDist[left] = temp;

   if ((left - begin) > 1)
   {
      qsortFitness(f, begin, left);
   }

   if ((lastPart - right) > 1)
   {
      qsortFitness(f, right, lastPart);
   }
}/* end qsortFitness() */

/******************************************************************************
compute_distance()

Compute the crowding distance of each particle in noDomF
******************************************************************************/
void MOPSOCD::compute_distance(unsigned int f)
{
   int i, max;

   max = 1;
   for (i = 1; i < m_NumNonDom - 1; i++)
   {
      m_CrowdDist[i] = m_CrowdDist[i] + (m_ArchiveFit[i + 1][f] - m_ArchiveFit[i - 1][f]);

      if (m_CrowdDist[max] < m_CrowdDist[i])
      {
         max = i;
      }
   }

   /* give maximum crowding distance value to the boundary points so that they are always selected */
   m_CrowdDist[0] = m_CrowdDist[0] + m_CrowdDist[max];
   m_CrowdDist[m_NumNonDom - 1] = m_CrowdDist[m_NumNonDom - 1] + m_CrowdDist[max];

}/* end compute_distance() */

/******************************************************************************
qsortCrowd()

Sort crowding distance values.
******************************************************************************/
void MOPSOCD::qsortCrowd(unsigned int begin, unsigned int lastPart)
{
   unsigned int left = begin + 1;
   unsigned int right = lastPart;
   double pivot = m_CrowdDist[begin];

   int k;
   double temp;
   double * tempP, * tempF;

   tempP = new double[m_NumVar];
   tempF = new double[m_NumFun];

   while (left < right)
   {

      if (m_CrowdDist[left] >= pivot) 
      {
         left++;
      }
      else
      {
         right--;

         /* exchange their crowding distance values */
         temp = m_CrowdDist[left];
         m_CrowdDist[left] = m_CrowdDist[right];
         m_CrowdDist[right] = temp;

         /* Exchange fitness positions of two particles in the array archiveFit */
         for (k = 0; k < m_NumFun; k++)
         {
            tempF[k] = m_ArchiveFit[left][k];
            m_ArchiveFit[left][k] = m_ArchiveFit[right][k];
         }

         for (k = 0; k < m_NumFun; k++)
         {
            m_ArchiveFit[right][k] = tempF[k];
         }

         /* Also exchange particle positions in archiveVar */
         for (k = 0; k < m_NumVar; k++)
         {
            tempP[k] = m_ArchiveVar[left][k];
            m_ArchiveVar[left][k] = m_ArchiveVar[right][k];
         }

         for (k = 0; k < m_NumVar; k++)
         {
            m_ArchiveVar[right][k] = tempP[k];
         }
      }/* end else() */
   }/* end while() */

   left--;

   /* Exchange fitness positions of two particles in the array archiveVar */
   for (k = 0; k < m_NumFun; k++)
   {
      tempF[k] = m_ArchiveFit[begin][k];
      m_ArchiveFit[begin][k] = m_ArchiveFit[left][k];
   }

   for (k = 0; k < m_NumFun; k++)
   {
      m_ArchiveFit[left][k] = tempF[k];
   }

   /* Also exchange particle positions in archiveVar */
   for (k = 0; k < m_NumVar; k++)
   {
      tempP[k] = m_ArchiveVar[begin][k];
      m_ArchiveVar[begin][k] = m_ArchiveVar[left][k];
   }

   for (k = 0; k < m_NumVar; k++)
   {
      m_ArchiveVar[left][k] = tempP[k];
   }

   /* Also exchange their crowding distance */
   temp = m_CrowdDist[begin];
   m_CrowdDist[begin] = m_CrowdDist[left];
   m_CrowdDist[left] = temp;

   if (left - begin > 1)
   {
      qsortCrowd(begin, left);
   }

   if (lastPart - right > 1)
   {
      qsortCrowd(right, lastPart);
   }

   delete [] tempP;
   delete [] tempF;
}/* end qsortCrowd() */

/******************************************************************************
compute_velocity()

Compute new velocity of each particle in the population
******************************************************************************/
void MOPSOCD::compute_velocity(void) 
{
   unsigned int top, gBest;
   int i, j;

   top = (unsigned int)((m_NumNonDom - 1) * 0.10);

   for (i = 0; i < m_PopSize; i++) 
   {
      gBest = MOPSOCD_RandomInt(0, top);

      for (j = 0; j < m_NumVar; j++)
      {
         m_Velocity[i][j] = 0.4 * m_Velocity[i][j] + /* W  * Vi */
                            1.0 * MOPSOCD_RandomDouble(0.0, 1.0) * (m_PbestsVar[i][j] - m_PopVar[i][j]) +     /* C1  * RandomDouble(0.0, 1.0) * (pBest - Xi) */
                            1.0 * MOPSOCD_RandomDouble(0.0, 1.0) * (m_ArchiveVar[gBest][j] - m_PopVar[i][j]); /* C2  * RandomDouble(0.0, 1.0) * (gBest - Xi) */          
      }/* end for() */
   }/* end for() */

   /* Calculate new positions of particles */
   for (i = 0; i < m_PopSize; i++)
   {
      for (j = 0; j < m_NumVar; j++)
      {
         m_PopVar[i][j] = m_PopVar[i][j] + m_Velocity[i][j];
      }
   }
}/* end compute_velocity() */

/******************************************************************************
mutate()

Mutation of particles adapted from MOPSO
******************************************************************************/
void MOPSOCD::mutate(unsigned int t)
{
   int i;
   int dimension = 0;
   double * minvalue, * maxvalue;
   double minvaluetemp, maxvaluetemp, range;
   double valtemp = 0;

   minvalue = new double[m_NumVar];
   maxvalue = new double[m_NumVar];

   get_ranges(minvalue, maxvalue);

   for (i = 0; i < m_PopSize; i++)
   {
      if (MOPSOCD_flip(pow(1.0 - (double)t / (m_MaxGen * m_ProbMut), 1.5)))
      {
         dimension = MOPSOCD_RandomInt(0, m_NumVar - 1);

         range = (maxvalue[dimension] - minvalue[dimension]) * pow(1.0 - (double)t / (m_MaxGen * m_ProbMut), 1.5) / 2;

         valtemp = MOPSOCD_RandomDouble(range, -range);

         if ((m_PopVar[i][dimension] - range) < minvalue[dimension])
         {
            minvaluetemp = minvalue[dimension];
         }
         else
         {
            minvaluetemp = m_PopVar[i][dimension] - range;
         }

         if (m_PopVar[i][dimension] + range > maxvalue[dimension])
         {
            maxvaluetemp = maxvalue[dimension];
         }
         else
         {
            maxvaluetemp = m_PopVar[i][dimension] + range;
         }

         m_PopVar[i][dimension] = MOPSOCD_RandomDouble(minvaluetemp, maxvaluetemp);
      }/* end if() */
   }/* end for() */

   delete [] minvalue;
   delete [] maxvalue;
}/* end mutate() */

/******************************************************************************
update_archive()

Insert new nondominated particles in pop into archive.
******************************************************************************/
void MOPSOCD::update_archive(void)
{
   int i, j, k;
   unsigned int bottom;

   /* for each particle in the population */
   for (k = 0; k < m_PopSize; k++) 
   {
      /* if particle in pop is nondominated */
      if (check_nondom(k) == 1) 
      {
         /* if memory is not yet full, insert particle */
         if (m_NumNonDom < m_ArchiveSize) 
         {
            i = m_NumNonDom;

            for (j = 0; j < m_NumVar; j++)
            {
               m_ArchiveVar[i][j] = m_PopVar[k][j];
            }

            for (j = 0; j < m_NumFun; j++)
            {
               m_ArchiveFit[i][j] = m_PopFit[k][j];
            }

            m_NumNonDom += 1;
         }
         else /* if memory is full, select particle to replace */
         {      
            /* Crowding distance computation and sorting */
            crowding();
            bottom = (unsigned int)((m_NumNonDom - 1)*.90);
            i = MOPSOCD_RandomInt(bottom, m_NumNonDom - 1);

            /* insert new particle to noDom */
            for (j = 0; j < m_NumVar; j++)
            {
               m_ArchiveVar[i][j] = m_PopVar[k][j];
            }

            for (j = 0; j < m_NumFun; j++)
            {
               m_ArchiveFit[i][j] = m_PopFit[k][j];
            }
         }/* end else() */
      }/* end if() */
   }/* end for() */
}/* end update_archive() */

/******************************************************************************
check_nondom()

Check for feasibility and nondomination of particles in pop and archive
******************************************************************************/
unsigned int MOPSOCD::check_nondom(unsigned int i)
{
   unsigned int sum;
   int h = 0;
   int j;
   double * archiveCons, * popCons;

   archiveCons = new double[m_NumVar];
   popCons = new double[m_NumVar];

   do 
   {
      sum = 0;
      for (j = 0; j < m_NumVar; j++)
      {
         archiveCons[j] = m_ArchiveVar[h][j];
         popCons[j] = m_PopVar[i][j];
      }

      /* if particle in archive has lesser contraint violations */
      if (check_constraints(archiveCons) < check_constraints(popCons))
      {
         delete [] archiveCons;
         delete [] popCons;
         return 0;
      }

      /* if particle in archive has more contraint violations, delete it */
      if (check_constraints(archiveCons) > check_constraints(popCons))
      {
         for (j = 0; j < m_NumVar; j++)
         {
            m_ArchiveVar[h][j] = m_ArchiveVar[m_NumNonDom - 1][j];
         }

         for (j = 0; j < m_NumFun; j++)
         {
            m_ArchiveFit[h][j] = m_ArchiveFit[m_NumNonDom - 1][j];
         }

         m_NumNonDom -= 1;
      }
      else
      {
         for (j = 0; j < m_NumFun; j++)
         {
            if (((m_ArchiveFit[h][j] < m_PopFit[i][j]) && (MOPSOCD_OPTIMIZATION_TYPE == 0)) ||
                ((m_ArchiveFit[h][j] > m_PopFit[i][j]) && (MOPSOCD_OPTIMIZATION_TYPE == 1)))
            {
               sum += 1;
            }/* end if() */
         }/* end for() */

         if (sum == m_NumFun)	/* If particle in archive dominates */
         {
            delete[] archiveCons;
            delete[] popCons;
            return 0;
         }
         else if (sum == 0) /* If particle in archive is dominated, delete it */
         {	
            for (j = 0; j < m_NumVar; j++)
            {
               m_ArchiveVar[h][j] = m_ArchiveVar[m_NumNonDom - 1][j];
            }

            for (j = 0; j < m_NumFun; j++)
            {
               m_ArchiveFit[h][j] = m_ArchiveFit[m_NumNonDom - 1][j];
            }
            m_NumNonDom -= 1;
         }
         else 
         {
            h += 1;
         }
      }/* end else() */
   } while (h < m_NumNonDom);

   delete[] archiveCons;
   delete[] popCons;
   return 1;
}/* end check_nondom() */

/******************************************************************************
get_ranges()

Get range values of variables
******************************************************************************/
void MOPSOCD::get_ranges(double *minvalue, double *maxvalue)
{
   int i;

   for (i = 0; i < m_NumVar; i++)
   {
      minvalue[i] = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLwrBnd();
      maxvalue[i] = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUprBnd();
   }
}/* end get_ranges() */

/******************************************************************************
maintain_particles()

Maintain particles in the population within the search space
******************************************************************************/
void MOPSOCD::maintain_particles(void)
{
   int i, j;
   double minvalue, maxvalue;

   for (i = 0; i < m_PopSize; i++) 
   {
      for (j = 0; j < m_NumVar; j++) 
      {
         minvalue = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetLwrBnd();
         maxvalue = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetUprBnd();

         /* If particle goes beyond minimum range value */
         if (m_PopVar[i][j] < minvalue)
         {

            /* Set it to minimum range value */
            m_PopVar[i][j] = minvalue;

            /* Change to opposite direction */
            m_Velocity[i][j] = -m_Velocity[i][j];

         }

         /* If particle goes beyond maximum range value */
         if (m_PopVar[i][j] > maxvalue)
         {
            /* Set it to maximum range value */
            m_PopVar[i][j] = maxvalue;

            /* Change to opposite direction */
            m_Velocity[i][j] = -m_Velocity[i][j];
         }
      }/* end for() */
   }/* end for() */
}/* maintain_particles() */

/******************************************************************************
update_pbests()

Update personal bests of particles in the population
******************************************************************************/
void MOPSOCD::update_pbests(void) 
{
   unsigned int sum, better;
   int i, j;

   for (i = 0; i < m_PopSize; i++) 
   {
      sum = 0;

      for (j = 0; j < m_NumFun; j++)
      {
         if (((m_PopFit[i][j] < m_PbestsFit[i][j]) && (MOPSOCD_OPTIMIZATION_TYPE == 0)) ||
             ((m_PopFit[i][j] > m_PbestsFit[i][j]) && (MOPSOCD_OPTIMIZATION_TYPE == 1)))
         {
            sum += 1;
         }
      }

      if (sum == m_NumFun) 
      {
         better = 0;
      }
      else 
      {
         if (sum == 0)
         {
            better = 1;
         }
         else
         {
            better = MOPSOCD_RandomInt(0, 1);
         }
      }

      if (better == 0)
      {
         for (j = 0; j < m_NumFun; j++)
         {
            m_PbestsFit[i][j] = m_PopFit[i][j];
         }
         for (j = 0; j < m_NumVar; j++)
         {
            m_PbestsVar[i][j] = m_PopVar[i][j];
         }
      }
   }/* end for() */
}/* end update_pbests() */

/******************************************************************************
save_results()

Write results to file
******************************************************************************/
void MOPSOCD::save_results(char *archiveName)
{
   int i, j;
   FILE *fp;

   /* Open file for writing */
   fp = fopen(archiveName, "w");

   for (i = 0; i < m_NumNonDom; i++) 
   {
      for (j = 0; j < m_NumFun; j++)
      {
         fprintf(fp, "%.6f ", m_ArchiveFit[i][j]);
      }
      fprintf(fp, "\n");
   }

   fclose(fp);
} /* end save_results() */

/******************************************************************************
MOPSOCD_Program()

Calibrate or optimize the model using MOPSOCD.
******************************************************************************/
void MOPSOCD_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("MOPSOCD", 1);
   MOPSOCD * TheAlg = new MOPSOCD(model);
   MEM_CHECK(TheAlg);

   if (model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end MOPSOCD_Program() */

