/******************************************************************************
File     : ParticleSwarm.h
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Particle Swarm Optimization (PSO) applies concepts of social behavior to
solve optimization problems. The PSO algorithm starts with a 'swarm' of 
particles (solutions) and "flies" this population through the design space in
search of the optimal solution. At each iteration, a given particle uses it's
own prior best solution (cognitive behavior) along with the current best 
solution of all particles (social behavior) to decide where to go next.

Version History
02-25-04    lsm   added copyright information and initial comments.
03-25-04    lsm   added PSO-LevMar hybrid
07-08-04    lsm   added parallel support, switched over to ParameterABC
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
11-18-04    lsm   Added convergence criteria, based on median fitness of swarm.
12-02-04    lsm   Added support for the seeding of particle swarm
10-19-05    lsm   Added support for linearly reducing the inertia weight to zero.
                  Replaced calls to rand() with MyRand(). Added support for status
                  file (used in grid computing).
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel. 
01-01-07    lsm   Latin Hypercube Sampling has been added as an option for 
                  initializing the PSO swarm. User's select this option by
                  including the following line in the PSO section of the 
                  input file:
                     InitPopulationMethod LHS
07-18-07    lsm   Added support for SuperMUSE
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "ParticleSwarm.h"
#include "LevenbergAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "StatsClass.h"
#include "LatinHypercube.h"
#include "QuadTree.h"
#include "SuperMUSE.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"
#include "SuperMuseUtility.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void ParticleSwarm::WarmStart(void)
{
   int i;
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);

   for(i = 0; i < np; i++)
   {
      m_pSwarm[0].x[i] = pbest[i];
      m_pSwarm[0].b[i] = pbest[i];
   }
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
ParticleSwarm::ParticleSwarm(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pSwarm = NULL;
   m_pStats = NULL;
   m_pTrees = NULL;
   m_pInit = NULL;
   m_TreeSize = 0;
   m_SwarmSize = 0;
   m_BestIdx = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumInit = 0;

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   m_pBuf    = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the PSO and it's member variables.
******************************************************************************/
void ParticleSwarm::Destroy(void)
{
   int i;
   for(i = 0; i < m_SwarmSize; i++)
   {
      delete [] m_pSwarm[i].v;
      delete [] m_pSwarm[i].x;
      delete [] m_pSwarm[i].b;
      delete [] m_pSwarm[i].cb;
      delete [] m_pSwarm[i].cx;
   }
   delete [] m_pSwarm;

   for(i = 0; i < m_NumInit; i++)
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;
   
   delete [] m_Fmedian;
   delete m_pStats;
   delete [] m_pTrees;
   m_SwarmSize = 0;
   m_BestIdx = 0;

   delete [] m_pMyBuf;
   delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   delete [] m_pBuf;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
CalcPSOMedian()

Calculates and returns the median objective function of the swarm. This 
parameter is used in the termination criteria of the PSO algorithm.
******************************************************************************/  
double ParticleSwarm::CalcPSOMedian(void)
{
   double med;
   int i;
    
   for(i = 0; i < m_SwarmSize; i++) { m_Fmedian[i] = m_pSwarm[i].fx; }
   med = CalcMedian(m_Fmedian, m_SwarmSize);

   return med;
}/* end CalcPSOMEdian() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PSO.
******************************************************************************/
void ParticleSwarm::Calibrate(void)
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

Minimize the objective function using PSO.
******************************************************************************/
void ParticleSwarm::Optimize(void)
{
   // char msg[DEF_STR_SZ];
   StatusStruct pStatus;
   int num, g;
   int id, lvl, idx;
   int i, j;
   double upr, lwr, r, range, sgn;
   double avg, x, pl, pg, r1, r2, v, vmin, median;
   double rval, init; //initial inertia
   double * pVals;
   LatinHypercube * pLHS = NULL;
   ParameterGroup * pGroup;
   
   InitFromFile(GetInFileName());

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if(id == 0)
   {
      WriteSetup(m_pModel, "Particle Swarm Optimization");
      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Convergence Value");
   }/* end if() */

   //allocate swarm
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();
   int nSpecial = pGroup->GetNumSpecialParams();

   NEW_PRINT("ParticleStruct", m_SwarmSize);
   m_pSwarm = new ParticleStruct[m_SwarmSize];
   MEM_CHECK(m_pSwarm);

   NEW_PRINT("double", m_SwarmSize);
   m_Fmedian = new double[m_SwarmSize];
   MEM_CHECK(m_Fmedian);

   for(i = 0; i < m_SwarmSize; i++)
   {
      NEW_PRINT("double", num);
      m_pSwarm[i].x = new double[num];

      NEW_PRINT("double", num);
      m_pSwarm[i].v = new double[num];

      NEW_PRINT("double", num);
      m_pSwarm[i].b = new double[num];

      NEW_PRINT("double", nSpecial);
      m_pSwarm[i].cb = new double[nSpecial];

      NEW_PRINT("double", nSpecial);
      m_pSwarm[i].cx = new double[nSpecial];

      m_pSwarm[i].n = num;
   }/* end for() */
   MEM_CHECK(m_pSwarm[i-1].b);

   //initialize swarm
   if(m_InitType == LHS_INIT)
   {
      NEW_PRINT("LatinHypercube", 1);
      pLHS = new LatinHypercube(num, m_SwarmSize);
      MEM_CHECK(pLHS);

      for(j = 0; j < num; j++)
      { 
         lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
         upr = pGroup->GetParamPtr(j)->GetUprBnd();
         pLHS->InitRow(j, lwr, upr);
      }/* end for() */
   }/* end if() */
 
   lvl = idx = 0;
   for(i = 0; i < m_SwarmSize; i++) //for each particle
   {
      //initial velocity is 0.00
      for(j = 0; j < num; j++){ m_pSwarm[i].v[j] = 0.00;}

      if(m_InitType == RANDOM_INIT)
      {
         for(j = 0; j < num; j++) //for each parameter
         {
            //generate a random between lower and upper bound
            lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
            upr = pGroup->GetParamPtr(j)->GetUprBnd();
            range = upr - lwr;
            r = (double)MyRand() / (double)MY_RAND_MAX;
            rval = (r * range) + lwr;
            m_pSwarm[i].x[j] = rval;
            m_pSwarm[i].b[j] = rval;         
         }/* end for() */
      }/* end if() */
      else if(m_InitType == QUAD_TREE_INIT)
      {
         //initialize quad trees if needed
         if(m_pTrees == NULL)
         {
            m_TreeSize = num;
            NEW_PRINT("QuadTree", m_TreeSize);
            m_pTrees = new QuadTree[m_TreeSize];
            for(j = 0; j < m_TreeSize; j++)
            { 
               lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
               upr = pGroup->GetParamPtr(j)->GetUprBnd();
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
         for(j = 0; j < num; j++)
         {
            m_pSwarm[i].x[j] = pVals[j];
            m_pSwarm[i].b[j] = pVals[j];
         }/* end for() */
         delete [] pVals;
      }/* end else if(QUAD_TREE_INIT) */
      else //LHS_INIT
      {
         for(j = 0; j < num; j++)
         { 
            rval = pLHS->SampleRow(j);
            m_pSwarm[i].x[j] = rval;
            m_pSwarm[i].b[j] = rval;
         }/* end for() */
      }/* end else() */
   }/* end for() */

   //seed swarm with pre-specified values
   for(i = 0; (i < m_NumInit) && (i < m_SwarmSize); i++)
   {
      for(j = 0; j < num; j++)
      {         
         m_pSwarm[i].x[j] = m_pInit[i][j];
         m_pSwarm[i].b[j] = m_pInit[i][j];
      }/* end for() */
   }/* end for() */

   //insert warm start solution, if desired
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }
   //insert extracted solution, if desired
   if(pGroup->CheckExtraction() == true)
   {
      pGroup->ReadParams(m_pSwarm[0].x);
      pGroup->ReadParams(m_pSwarm[0].b);
   }/* end if() */

   delete pLHS;

   //evaluate swarm, possibly in parallel
   m_CurGen = 0;
   EvaluateSwarm();

   //perform intermediate bookkeeping
   m_pModel->Bookkeep(false);

   for(i = 0; i < m_SwarmSize; i++)
   {
	   m_pSwarm[i].fb = m_pSwarm[i].fx;
      for(int iSpecial = 0; iSpecial < nSpecial; iSpecial++)
      {
	      m_pSwarm[i].cb[iSpecial] = m_pSwarm[i].cx[iSpecial];
      }
      //sync best with current, in case parameter corrections were made
      for(j = 0; j < num; j++)
      {
         m_pSwarm[i].b[j] =  m_pSwarm[i].x[j];
      }
   }

   /* --------------------------------------------
   enable special parameters now that local best 
   is initialized for each particle
   -------------------------------------------- */
   pGroup->EnableSpecialParams();

   //determine the best particle and average value
   m_BestIdx = 0;
   avg = 0.00;
   m_Best = m_pSwarm[0].fb;
   for(i = 0; i < m_SwarmSize; i++)
   {
      avg += m_pSwarm[i].fx;
      if(m_Best > m_pSwarm[i].fx)
      {
         m_Best = m_pSwarm[i].fx;
         m_BestIdx = i;
      }/* end if() */
   }/* end if() */
   avg /= m_SwarmSize;
   median = CalcPSOMedian();
   //current convergence value
   m_CurStop = fabs((median - m_Best)/median);

   if(id == 0)
   {
      //write initial config.
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
      WriteRecord(m_pModel, 0, m_Best, m_CurStop);
      pStatus.curIter = 0;
      pStatus.maxIter = m_MaxGens;
      pStatus.pct = 0.00;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
   }/* end if() */

   init = m_Inertia;
   //main optimization loop   
   for(g = 0; g < m_MaxGens; g++)
   {
      pStatus.curIter = m_CurGen = g+1;
      if(IsQuit() == true){ break;}
      if(m_CurStop < m_StopVal){ pStatus.pct = 100.00; break;}

      if(id == 0)
      {
         //update velocities, parameters and objective functions
         for(i = 0; i < m_SwarmSize; i++) //for each particle
         {
            for(j = 0; j < num; j++) //for each parameter
            {
               //intermediary variables
               x = m_pSwarm[i].x[j];
               pl = m_pSwarm[i].b[j];
               pg = m_pSwarm[m_BestIdx].b[j];

               //random weights
               r1 = (double)MyRand() / (double)MY_RAND_MAX;
               r2 = (double)MyRand() / (double)MY_RAND_MAX;

               //revised velocity
               v = m_pSwarm[i].v[j];
               v = m_Constrict*((m_Inertia*v) + m_c1*r1*(pl-x) + m_c2*r2*(pg-x));

   			   //assign minimum perturbation to prevent stagnation
	   		   if(strcmp(pGroup->GetParamPtr(j)->GetType(), "real") == 0)
		   		 vmin = (0.01*fabs(x))/(g+1); 
			      else
				    vmin = 0.50;

               if(fabs(v) < vmin)
               {
                  //adjust randomized minimum velocity
                  sgn = (double)MyRand() / (double)MY_RAND_MAX; //random direction
                  if(sgn >= 0.50) v = +((1.00+r1)*vmin);
                  else            v = -((1.00+r2)*vmin);
               }
               m_pSwarm[i].v[j] = v;
            
               //revised position
               m_pSwarm[i].x[j] = x + v;
            }/* end for() */

            /*-----------------------------------------------------------
            Constrain revised position to stay within parameter limits
            but be sure to preserve angle (i.e. direction of movement)
            -----------------------------------------------------------*/
            double dx_min, dx_old, dx_new, dx_frac;
            dx_min = 1.00;
            for(j = 0; j < num; j++) //for each parameter
            {
               lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
               upr = pGroup->GetParamPtr(j)->GetUprBnd();

               v = m_pSwarm[i].v[j];
               x = m_pSwarm[i].x[j] - v; //original position
               //compute fractional move, store most restrictive fraction
               if(m_pSwarm[i].x[j] > upr)
               {
                  dx_old = v;
                  dx_new = 0.5*(upr-x);
                  dx_frac = fabs(dx_new/dx_old); //relative change
                  if(dx_frac < dx_min) dx_min = dx_frac;
                  m_NumUprViols++;
               }
               if(m_pSwarm[i].x[j] < lwr)
               {
                  dx_old = v;
                  dx_new = 0.5*(lwr-x);
                  dx_frac = fabs(dx_new/dx_old); //relative change
                  if(dx_frac < dx_min) dx_min = dx_frac;
                  m_NumLwrViols++;
               }
            }/* end for() */
            for(j = 0; j < num; j++) //for each parameter
            {
               v = m_pSwarm[i].v[j];
               x = m_pSwarm[i].x[j] - v; //original position      
               m_pSwarm[i].v[j] *= dx_min; //revised velocity
               m_pSwarm[i].x[j] = x + (v*dx_min); //revised position
            }
         }/* end for() */
      } /* end if() */

      //evaluate swarm, possibly in parallel
      EvaluateSwarm();

      //reduce inertia
      if(m_LinRedFlag == true) //linearly reducing to zero
      {
         m_Inertia = init;
         m_RedRate = (double)g/(double)m_MaxGens;
      }
      m_Inertia *= (1.00 - m_RedRate);

      //FILE * pWgt = fopen("PSO.txt", "a");
      //fprintf(pWgt, "Inertia = %E, Cognitive = %E, Social = %E\n", m_Inertia, m_c1, m_c2);
      //fclose(pWgt);
   
      //revise avg, and local and global best
      avg = 0.00;
      for(i = 0; i < m_SwarmSize; i++)
      {
         avg += m_pSwarm[i].fx;
         //revise local best
         if(m_pSwarm[i].fx < m_pSwarm[i].fb)
         {
            for(j = 0; j < num; j++){ m_pSwarm[i].b[j] = m_pSwarm[i].x[j];}
            m_pSwarm[i].fb = m_pSwarm[i].fx;

            for(int iSpecial = 0; iSpecial < nSpecial; iSpecial++)
            {
			      m_pSwarm[i].cb[iSpecial] = m_pSwarm[i].cx[iSpecial];
            }
         }/* end if() */
         //revise global best
         if(m_Best > m_pSwarm[i].fx)
         {
            m_Best = m_pSwarm[i].fx;
            m_BestIdx = i;
         }/* end if() */
      }/* end for() */

      avg /= m_SwarmSize;
      median = CalcPSOMedian();
      //current convergence value
      m_CurStop = fabs((median - m_Best)/median);
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);

      if(id == 0){ WriteRecord(m_pModel, (g+1), m_Best, m_CurStop);}
      pStatus.pct = ((float)100.00*(float)(g+1))/(float)m_MaxGens;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */

   m_Inertia = init; //reset inertia

   //place model at optimal prameter set
   pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
   m_pModel->Execute();

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   { 
      WriteOptimal(m_pModel, m_Best);
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
void ParticleSwarm::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Particle Swarm Optimization\n");
   fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
   fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
   fprintf(pFile, "Max Generations         : %d\n", m_MaxGens);
   fprintf(pFile, "Actual Generations      : %d\n", m_CurGen);
   fprintf(pFile, "Swarm Size              : %d\n", m_SwarmSize);
   fprintf(pFile, "Constriction Factor     : %.2lf\n", m_Constrict);  
   fprintf(pFile, "Cognitive Weight        : %.2lf\n", m_c1);
   fprintf(pFile, "Social Weight           : %.2lf\n", m_c2);
   fprintf(pFile, "Inertia Weight          : %.2lf\n", m_Inertia);
   
   fprintf(pFile, "Inertia Reduction Rate  : ");
   if(m_LinRedFlag == true) fprintf(pFile, "Linear reduction to zero\n");
   else                     fprintf(pFile, "%.2lf\n", m_RedRate);

   fprintf(pFile, "Initialization Method   : ");
   if(m_InitType == RANDOM_INIT){ fprintf(pFile, "Random\n");}
   else if(m_InitType == QUAD_TREE_INIT){ fprintf(pFile, "Quad-Tree\n");}
   else if(m_InitType == LHS_INIT){ fprintf(pFile, "Latin Hypercube Sampling\n");}
   else { fprintf(pFile, "Unknown\n");}

   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());      
   fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

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
EvaluateSwarm()

Evaluates the objective function of each particle in the swarm.
******************************************************************************/
void ParticleSwarm::EvaluateSwarm(void)
{
   static double a = 0.00;
   int i, n, id;   
   ParameterGroup * pGroup;
   double val;

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   
   if(n == 1) //serial execution
   {
      if(IsSuperMUSE() == false)
      {
         WriteInnerEval(WRITE_PSO, m_SwarmSize, '.');
         pGroup = m_pModel->GetParamGroupPtr();
         for(i = 0; i < m_SwarmSize; i++) 
         { 
            WriteInnerEval(i+1, m_SwarmSize, '.');
            MakeParameterCorrections(m_pSwarm[i].x, m_pSwarm[m_BestIdx].b, m_pModel->GetParamGroupPtr()->GetNumParams(), a);

			   //let special parameters know about local best
			   pGroup->ConfigureSpecialParams(m_pSwarm[i].fb, m_pSwarm[i].cb);

            val = m_pModel->Execute();
            a += 1.00/(double)(m_SwarmSize*(m_MaxGens+1));
            m_pSwarm[i].fx = val;
			   pGroup->GetSpecialConstraints(m_pSwarm[i].cx);
         }
         WriteInnerEval(WRITE_ENDED, m_SwarmSize, '.');
      }
      else
      {
         EvalSwarmSuperMUSE();
      }
   }/* end if() */
   else /* parallel execution */
   {
      if(id == 0)
      {
         for(i = 0; i < m_SwarmSize; i++) 
         { 
            MakeParameterCorrections(m_pSwarm[i].x, m_pSwarm[m_BestIdx].b, m_pModel->GetParamGroupPtr()->GetNumParams(), a);
            a += 1.00/(double)(m_SwarmSize*(m_MaxGens+1));
         }/* end for() */
      }/* end if() */

      BcastSwarm();
      EvalSwarmParallel();      
   }/* end else() */
} /* end EvaluateSwarm() */

/*************************************************************************************
MakeParameerCorrections()
*************************************************************************************/
void ParticleSwarm::MakeParameterCorrections(double * x, double * xb, int n, double a)
{
   double lwr, upr;
   ParameterGroup * pParamGroup;
	pParamGroup = m_pModel->GetParamGroupPtr(); 

   for(int k = 0; k < n; k++)
   {
      lwr=m_pModel->GetParamGroupPtr()->GetParamPtr(k)->GetLwrBnd();
      upr=m_pModel->GetParamGroupPtr()->GetParamPtr(k)->GetUprBnd();
      x[k]=TelescopicCorrection(lwr, upr, xb[k], a, x[k]);
   }
   pParamGroup->WriteParams(x); 		

   //inerface with expert judgement module
   m_pModel->PerformParameterCorrections();
   for(int i = 0; i < n; i++)
   {
      x[i] = pParamGroup->GetParamPtr(i)->GetEstVal();
   }/* end for() */
}/* enad MakeParameterCorrections() */

/******************************************************************************
BcastSwarm()

When in parallel, only the master computes the swarm movement. All the other 
processors just compute the objeective functions. The BcastSwarm() routine is 
called upon to broadcast  the current particle swarm from the master processor 
to all of the slave processors.
******************************************************************************/
void ParticleSwarm::BcastSwarm(void)
{
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix
   pop_size = m_SwarmSize;
   num_vars = m_pSwarm[0].n;

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
         m_pBuf[idx] = m_pSwarm[j].x[i];
      }/* end for() */
   }/* end for() */

   //broadcast the flattened matrix
   MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //use the flattened matrix to fill swarm
   for(i = 0; i < num_vars; i++)
   {
      for(j = 0; j < pop_size; j++)
      {
         
         idx = (num_vars)*j + i;
         m_pSwarm[j].x[i] = m_pBuf[idx];
      }/* end for() */
   }/* end for() */
}/* end BcastSwarm() */

/******************************************************************************
EvalSwarmParallel()

Compute objective function of entire particle swarm in parallel. Each processor 
evaluates a predetermined number of swarm particles, based on their processor id.
******************************************************************************/
void ParticleSwarm::EvalSwarmParallel(void)
{    
   int i ,j, num_procs, id, bufsize, idx;
   ParameterGroup * pGroup;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
   bufsize = (m_SwarmSize/num_procs) + 1;

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];

      NEW_PRINT("double", bufsize);
      m_pTmpBuf = new double[bufsize];

      NEW_PRINT("double", m_SwarmSize);
      m_pBigBuf = new double[m_SwarmSize];
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_SwarmSize; i++) 
   { 
      if((i % num_procs) == id)
      {          
         pGroup->WriteParams(m_pSwarm[i].x);

		   //let special parameters know about local best
		   pGroup->ConfigureSpecialParams(m_pSwarm[i].fb, m_pSwarm[i].cb);

         m_pMyBuf[j] = m_pModel->Execute();
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

         if(idx < m_SwarmSize)
         {
            m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            m_pMyBuf[j] = m_pTmpBuf[j]; //restore myBuf...clobbered by bcast
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //stuff results into swarm
   for(i = 0; i < m_SwarmSize; i++)
   {
      m_pSwarm[i].fx = m_pBigBuf[i];
   }/* end for() */
}/* end EvalSwarmParallel() */

/******************************************************************************
EvalSwarmSuperMUSE()

Compute objective functions of the swarm using SuperMUSE. This routine 
interfaces with the RepeatTasker SuperMUSE program, which assigns model 
evaluations to SuperMUSE clients on a first-come-first-served basis.
******************************************************************************/
void ParticleSwarm::EvalSwarmSuperMUSE(void)
{  
   double val;
   bool bOk; 
   int i;
   ParameterGroup * pGroup;
   SuperMUSE * pSMUSE = GetSuperMusePtr();

   /* ----------------------------------------------------------------
   Generate task file that describes the desired parallel evaluations.
   This is somewhat analogous to the BcastGrid() operation used
   for MPI-parallel operations. Write the parameter values of each 
   swarm member as entries in the task file.

   Entries are first accumlated into a temp file to prevent the 
   SuperMUSE RepeatTasker program from prematurely processing the task
   file.
   ---------------------------------------------------------------- */
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_SwarmSize; i++)
   {
      //stuff the parameter group with values
      pGroup->WriteParams(m_pSwarm[i].x);
         
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
      EvaluateSwarm();
   }
   else //SuperMUSE was successful
   {
      for(i = 0; i < m_SwarmSize; i++)
      {
         /* -----------------------------------------------
         Stuff the parameter group with ith swarm 
         member. This ensures that each objective function 
         gets associated with the correct parameter values.
         ------------------------------------------------ */
         pGroup->WriteParams(m_pSwarm[i].x);

         //stuff i-th result into chromosome pool
         val = pSMUSE->GatherResult(i);
         m_pSwarm[i].fx = val;
		   pGroup->GetSpecialConstraints(m_pSwarm[i].cx);
      }/* end for() */
   }/* end else() */
}/* end EvalSwarmSuperMUSE() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void ParticleSwarm::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   int i, j, k, num;
   char * pTok;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];

   m_StopVal = 0.001;
   m_SwarmSize = 20;
   m_MaxGens = 50;
   m_Constrict = 1.00;
   m_c1 = 2.00;
   m_c2 = 2.00;
   m_Inertia = 1.2;   
   m_RedRate = 0.10;
   m_LinRedFlag = false;
   m_InitType = RANDOM_INIT;

   //read in PSO configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open PSO config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginParticleSwarm", pFileName) == true)
   {
      FindToken(pFile, "EndParticleSwarm", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginParticleSwarm", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndParticleSwarm") == NULL)
      {         
         if(strstr(line, "SwarmSize") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_SwarmSize); 
         }/*end else if() */         
         else if(strstr(line, "NumGenerations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxGens);
         }
         else if(strstr(line, "ConstrictionFactor") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Constrict);
         }
         else if(strstr(line, "CognitiveParam") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_c1);
         }
         else if(strstr(line, "SocialParam") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_c2);
         }
         else if(strstr(line, "InertiaWeight") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Inertia);
         }
         else if(strstr(line, "InertiaReductionRate") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            MyStrLwr(tmp2);
            if(strcmp(tmp2, "linear") == 0){ m_LinRedFlag = true;}
            else{sscanf(line, "%s %lf", tmp, &m_RedRate);}
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
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   /* initialize some or all swarm members to specied values */
   rewind(pFile);
   if(CheckToken(pFile, "BeginInitParams", pFileName) == true)
   {
      FindToken(pFile, "EndInitParams", pFileName);
      rewind(pFile);

      //allocate space for the parameter list
      num = m_pModel->GetParamGroupPtr()->GetNumParams();

      //count the number of entries
      FindToken(pFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      m_NumInit = 0;
      while(strstr(line, "EndInitParams") == NULL)
      {
         m_NumInit++;
         line = GetNxtDataLine(pFile, pFileName);
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
      FindToken(pFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      i = 0;
      while(strstr(line, "EndInitParams") == NULL)
      {
         pTok = line;
         //extract values, one-by-one, making any necessary conversions
         for(k = 0; k < num; k++)
         {
            j = ExtractString(pTok, tmp);
            j = ValidateExtraction(j, k, num, "PSO::InitFromFile()");
            pTok += j;            
            m_pInit[i][k] = m_pModel->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
         }/* end for() */                  
         i++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
   }/* end if() */

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
PSO_Program()

Calibrate or optimize the model using PSO.
******************************************************************************/
void PSO_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("ParticleSwarm", 1);
   ParticleSwarm * PSO = new ParticleSwarm(model);
   MEM_CHECK(PSO);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { PSO->Calibrate(); }
   else { PSO->Optimize(); }

   delete PSO;
   delete model;
} /* end PSO_Program() */

/******************************************************************************
PSO_LEVMAR_Program()

Calibrate the model using PSO-Levenberg-Marquardt hybrid.
******************************************************************************/
void PSO_LEVMAR_Program(int argC, StringType argV[])
{
   int id;
   char file1[DEF_STR_SZ], file2[DEF_STR_SZ];
   char * gmlStr = (char *)"_GML";
   char * psoStr = (char *)"_PSO";
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("ParticleSwarm", 1);
   ParticleSwarm * PSO = new ParticleSwarm(model);
   MEM_CHECK(PSO);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) 
   { 
      SetIterationResidualsPrefix(psoStr,0);

      PSO->Calibrate();
      delete PSO;

      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      sprintf(file1, "OstOutputPSO%d.txt", id);
      sprintf(file2, "OstOutput%d.txt", id);
      remove(file1);
      rename(file2, file1);

      sprintf(file1, "OstModelPSO%d.txt", id);
      sprintf(file2, "OstModel%d.txt", id);
      remove(file1);
      rename(file2, file1);

      sprintf(file1, "OstErrorsPSO%d.txt", id);
      sprintf(file2, "OstErrors%d.txt", id);
      remove(file1);
      rename(file2, file1);

      if(id == 0)
      {
         sprintf(file1, "OstStatusPSO%d.txt", id);
         sprintf(file2, "OstStatus%d.txt", id);
         remove(file1);
         rename(file2, file1);
      }

      NEW_PRINT("LevenbergAlgorithm", 1);
      LevenbergAlgorithm * LA = new LevenbergAlgorithm(model, false);
      MEM_CHECK(LA);
      SetIterationResidualsPrefix(gmlStr, 0);
      SetTrialNumber(1);
      LA->Calibrate(); 
      delete LA;
   }
   else 
   {
      printf("Hybrid GML-PSO algorithm can only be used for calibration.\n"); 
   }
   
   delete model;
}/* end PSO_LEVMAR_Program() */
