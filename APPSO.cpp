/******************************************************************************
File     : APPSO.cpp
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Asynchronous Parallel Particle Swarm Optimization (APPSO).

A parallel version of PSO based on the asynchronous master-slave approach
of the PDDS algorithm.

Version History
12-30-14    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "APPSO.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "LatinHypercube.h"
#include "StatsClass.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"

#define APPSO_DO_WORK   (101)
#define APPSO_STOP_WORK (102)

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void APPSO::WarmStart(void)
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
APPSO::APPSO(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pSwarm = NULL;
   m_pStats = NULL;
   m_pInit = NULL;
   m_Assignments = NULL;
   m_Fmedian = NULL;
   m_SwarmSize = 0;
   m_BestIdx = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumInit = 0;
   m_id = 0;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the PSO and it's member variables.
******************************************************************************/
void APPSO::Destroy(void)
{
   int i;
   if(m_id == 0)
   {
      for(i = 0; i < m_SwarmSize; i++)
      {
         delete [] m_pSwarm[i].v;
         delete [] m_pSwarm[i].x;
         delete [] m_pSwarm[i].b;
         delete [] m_pSwarm[i].cb;
         delete [] m_pSwarm[i].cx;
      }
   }
   delete [] m_pSwarm;

   for(i = 0; i < m_NumInit; i++)
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;
   
   delete [] m_Fmedian;
   delete [] m_Assignments;
   delete m_pStats;
   m_SwarmSize = 0;
   m_BestIdx = 0;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
CalcPSOMedian()

Calculates and returns the median objective function of the swarm. This 
parameter is used in the termination criteria of the PSO algorithm.
******************************************************************************/  
double APPSO::CalcPSOMedian(void)
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
void APPSO::Calibrate(void)
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
void APPSO::Optimize(void)
{
   StatusStruct pStatus;
   int num, g;
   int id;
   int i, j, nprocs, nSpecial;
   double upr, lwr, r, range, sgn;
   double avg, x, pl, pg, r1, r2, v, vmin, median;
   double rval, init; //initial inertia
   ParameterGroup * pGroup;

   InitFromFile(GetInFileName());

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   m_id = id;

   /* sanity check */
   if(nprocs < 2)
   {
      LogError(ERR_ABORT, "APPSO requires at least 2 processors");
      ExitProgram(0);  
   }

   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();
   nSpecial = pGroup->GetNumSpecialParams();

   m_Assignments = new int[nprocs];

   /* ===============================================================================
   Master processor allocates entire swarm and initializes swarm entries
   =============================================================================== */
   if(id == 0)
   {
      WriteSetup(m_pModel, "Particle Swarm Optimization");
      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Convergence Value");

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

      //random swarm initialization
      for(i = 0; i < m_SwarmSize; i++) //for each particle
      {
         //initial velocity is 0.00
         for(j = 0; j < num; j++){ m_pSwarm[i].v[j] = 0.00;}

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

      //read in best result from previous run, if desired
      if(m_pModel->CheckWarmStart() == true)
      {
         WarmStart();
      }
      // extract initial values, if desired
      if(pGroup->CheckExtraction() == true)
      {
         pGroup->ReadParams(m_pSwarm[0].x);
         pGroup->ReadParams(m_pSwarm[0].b);      
      }
   }/* end if(master) */
   /* ===============================================================================
   Slave processors allocates a single particle
   =============================================================================== */
   else
   {
      NEW_PRINT("ParticleStruct", 1);
      m_pSwarm = new ParticleStruct;
      MEM_CHECK(m_pSwarm);

      NEW_PRINT("double", num);
      m_pSwarm->x = new double[num];

      NEW_PRINT("double", num);
      m_pSwarm->v = new double[num];

      NEW_PRINT("double", num);
      m_pSwarm->b = new double[num];

      NEW_PRINT("double", nSpecial);
      m_pSwarm->cb = new double[nSpecial];

      NEW_PRINT("double", nSpecial);
      m_pSwarm->cx = new double[nSpecial];

      m_pSwarm->n = num;      
      MEM_CHECK(m_pSwarm->b);
   }/* end else() */

   //evaluate swarm, asynchronously and in parallel
   EvaluateSwarm(id, nprocs, 0);

   if(id == 0)
   {
      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);

      //assign initial best positions
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

      //write initial config.
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
      WriteRecord(m_pModel, 0, m_Best, median);
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
      } /* end if(master) */

      //evaluate swarm, possibly in parallel
      EvaluateSwarm(id, nprocs, g+1);

      if(id == 0)
      {
         //reduce inertia
         if(m_LinRedFlag == true) //linearly reducing to zero
         {
            m_Inertia = init;
            m_RedRate = (double)g/(double)m_MaxGens;
         }
         m_Inertia *= (1.00 - m_RedRate);
  
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

         pGroup->WriteParams(m_pSwarm[m_BestIdx].b);

         if(id == 0){ WriteRecord(m_pModel, (g+1), m_Best, median);}
         pStatus.pct = ((float)100.00*(float)(g+1))/(float)m_MaxGens;
         pStatus.numRuns = (g+1)*m_SwarmSize;
         WriteStatus(&pStatus);

         //perform intermediate bookkeeping
         m_pModel->Bookkeep(false);
      }/* end if(master) */
   }/* end for() */

   m_Inertia = init; //reset inertia

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   { 
      //place model at optimal prameter set
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
      m_pModel->Execute();

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
void APPSO::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Asynchronous Parallel Particle Swarm Optimization\n");
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

   fprintf(pFile, "Initialization Method   : Random\n");

   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());      
   fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
EvaluateSwarm()

Evaluates the objective function of each particle in the swarm.
******************************************************************************/
void APPSO::EvaluateSwarm(int id, int nprocs, int gen)
{
   MPI_Status mpi_status;
   static double a = 0.00;
   int i, ii, signal, num, sid, num_recv, nstops, nSpecial, nslaves, nxtsid;
   double * f;
   bool bDone = false;
   ParameterGroup * pGroup =  m_pModel->GetParamGroupPtr();
   bool bSynch = SynchReceives();

   num = pGroup->GetNumParams();
   nSpecial = pGroup->GetNumSpecialParams();
   f = new double[nSpecial+1];

   nstops = 0;
   nslaves = nprocs - 1;
   nxtsid = 0;

   //master "primes pump" with initial set of work
   if(id == 0)
   {
      //output banner
      WriteInnerEval(WRITE_PSO, m_SwarmSize, '.');

      //adjust parameters using meta heuristics and expert judgement
      for(i = 0; i < m_SwarmSize; i++) 
      { 
         MakeParameterCorrections(m_pSwarm[i].x, m_pSwarm[m_BestIdx].b, num, a);
         a += 1.00/(double)(m_SwarmSize*(m_MaxGens+1));
      }

      //assign initial work to slaves
      for(i = 1; i < nprocs; i++)
      {
         if(i <= m_SwarmSize)
         {
            m_Assignments[i] = i-1;
            // send work to slave
            signal = APPSO_DO_WORK;
            MPI_Send(&signal,1,MPI_INT,i,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
            MPI_Send(&(m_pSwarm[i-1].x[0]), num, MPI_DOUBLE, i, MPI_DATA_TAG, MPI_COMM_WORLD);
            MPI_Send(&(m_pSwarm[i-1].fb), 1, MPI_DOUBLE, i, MPI_DATA_TAG, MPI_COMM_WORLD);
            MPI_Send(&(m_pSwarm[i-1].cb[0]), nSpecial, MPI_DOUBLE, i, MPI_DATA_TAG, MPI_COMM_WORLD);
         }
         else
         {
            signal = APPSO_STOP_WORK;
            MPI_Send(&signal,1,MPI_INT,i,MPI_REQUEST_TAG,MPI_COMM_WORLD);
            nstops++;
         }
      }/* end for(each slave) */
   }/* end if(master) */

   num_recv = 0;
   while(bDone == false)
   {
      if(id == 0)
      {
         //receive result from slave and process
         if(bSynch == true)
         {
            sid = nxtsid + 1;
            nxtsid = (nxtsid + 1) % nslaves;
         }
         else
         {
            sid = MPI_ANY_SOURCE;
         }
         MPI_Recv(f, nSpecial+1, MPI_DOUBLE, sid, MPI_RESULTS_TAG, MPI_COMM_WORLD, &mpi_status);
         num_recv++;
         sid = mpi_status.MPI_SOURCE;

         /*
         FILE * pLog = fopen("OstMessages0.txt", "a");
         fprintf(pLog, "sid = %02d\n", sid);
         fclose(pLog);
         */

         WriteInnerEval(num_recv, m_SwarmSize, '.');
         ii = m_Assignments[sid];
         m_pSwarm[ii].fx = f[0];
         for(int iSpecial = 0; iSpecial < nSpecial; iSpecial++)
         {
   	      m_pSwarm[ii].cx[iSpecial] = f[1+iSpecial];
         }

         //assign more work
         if(i <= m_SwarmSize)
         {
            m_Assignments[sid] = i-1;
            // send work to slave
            signal = APPSO_DO_WORK;
            MPI_Send(&signal,1,MPI_INT,sid,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
            MPI_Send(&(m_pSwarm[i-1].x[0]), num, MPI_DOUBLE, sid, MPI_DATA_TAG, MPI_COMM_WORLD);
            MPI_Send(&(m_pSwarm[i-1].fb), 1, MPI_DOUBLE, sid, MPI_DATA_TAG, MPI_COMM_WORLD);
            MPI_Send(&(m_pSwarm[i-1].cb[0]), nSpecial, MPI_DOUBLE, sid, MPI_DATA_TAG, MPI_COMM_WORLD);
            i++;
         }
         else // send stop work message to the slave
         {
            signal = APPSO_STOP_WORK;
            MPI_Send(&signal,1,MPI_INT,sid,MPI_REQUEST_TAG,MPI_COMM_WORLD);
            nstops++;
            if(nstops == (nprocs-1)) //quit when all slaves have been told to stop
            {
               WriteInnerEval(WRITE_ENDED, m_SwarmSize, '.');
               bDone = true;
            }
         }
      }/* end if(master) */
      else /* slave processing */
      {
         MPI_Recv(&signal,1,MPI_INT,0,MPI_REQUEST_TAG,MPI_COMM_WORLD, &mpi_status); 
         if(signal == APPSO_DO_WORK)
         {
            num_recv++;
            MPI_Recv(&(m_pSwarm->x[0]), num, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Recv(&(m_pSwarm->fb), 1, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
            MPI_Recv(&(m_pSwarm->cb[0]), nSpecial, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
            if(num_recv == 1)
            {
               pGroup->EnableSpecialParams();
            }
            m_pModel->GetParamGroupPtr()->WriteParams(m_pSwarm->x);
            //let special parameters know about local best
	         pGroup->ConfigureSpecialParams(m_pSwarm->fb, m_pSwarm->cb);
            f[0] = m_pModel->Execute();
            pGroup->GetSpecialConstraints(&(f[1]));

            MPI_Send(f, 1+nSpecial, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);
         }/* end if() */
         else
         {
            bDone = true;
         }
      }/* end else(slave) */
   }/* end while() */

   //synch up processors
   MPI_Barrier(MPI_COMM_WORLD);

   delete [] f;
 } /* end EvaluateSwarm() */

/*************************************************************************************
MakeParameerCorrections()
*************************************************************************************/
void APPSO::MakeParameterCorrections(double * x, double * xb, int n, double a)
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
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void APPSO::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   int i, j, k, num;
   char * pTok;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];

   m_SwarmSize = 20;
   m_MaxGens = 50;
   m_Constrict = 1.00;
   m_c1 = 2.00;
   m_c2 = 2.00;
   m_Inertia = 1.2;   
   m_RedRate = 0.10;
   m_LinRedFlag = false;

   //read in PSO configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open APPSO config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginAPPSO", pFileName) == true)
   {
      FindToken(pFile, "EndAPPSO", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginAPPSO", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndAPPSO") == NULL)
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
            j = ValidateExtraction(j, k, num, "APPSO::InitFromFile()");
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
APPSO_Program()

Calibrate or optimize the model using PSO.
******************************************************************************/
void APPSO_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("APPSO", 1);
   APPSO * TheAlg = new APPSO(model);
   MEM_CHECK(TheAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end APPSO_Program() */

