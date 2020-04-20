/******************************************************************************
File      : SAAlgorithm.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

An implementation of the simulated annealing algorithm for continuously 
varying parameters.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   adjusted melt operation
                  outputs obj. func. str.
07-08-04    lsm   switched to ParameterABC
                  WriteSetup() is now the first action of algorithm
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
09-01-04    lsm   Algorithm now based on: Vanderbilt and Louie. 1984. "A Monte 
                  Carlo Simulated Annealing Approach to Optimization over 
                  Continuous Variables". Journal of Computational Physics. 
                  vol. 56, pg. 259-271.
11-18-04    lsm   Added convergence criteria, based on median F() of inner loop.
11-30-04    lsm   Replaced calls to time() with MyTime()
10-21-05    lsm   Switched to homegrown implementation that is easier to study
                  than the Vanderbilt and Louie implementation
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel. 
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "SAAlgorithm.h"
#include "Model.h"
#include "ModelBackup.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "TelescopingBounds.h"
#include "StatsClass.h"

#include "Utility.h"
#include "StatUtility.h"
#include "WriteUtility.h"
#include "Exception.h"

#define APSA_DO_WORK   (101)
#define APSA_STOP_WORK (102)

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void SAAlgorithm::WarmStart(void)
{
   int i;
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);

   /* -------------------------------------------------------------------------
   Depending on serial or parallel execution, the initial parameter set will 
   come from either m_pBest or the parameter group itself. So stuff the warm
   start info. into both locations.
   ----------------------------------------------------------------------- */
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   for(i = 0; i < np; i++)
   {
      m_pBest[i] = pbest[i];
   }
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Initializes parameters, reading user-specified input, if available.
******************************************************************************/
SAAlgorithm::SAAlgorithm(ModelABC * pModel)
{
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char tfnl[DEF_STR_SZ];
   char trans[DEF_STR_SZ];
   int numParams;
   ParameterGroup * pGroup;

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();
   m_pStats = NULL;
   m_pMelts = NULL;
   m_Finner = NULL;

   //init. everything to reasonable defaults
   m_InitProb = m_CurProb = -1.00;
   m_StopVal = 0.001;
   m_CurStop = 1.00;
   m_NumOuter = 0;
   m_MaxOuter = 20;
   m_MaxInner = 10;
   m_NumMelts = 100;
   m_InitTemp = m_CurTemp = m_FinalTemp = 10.00;
   m_FinalTempMethod = TMETHD_NORM;
   m_TransitionMethod = TRANS_GAUSS;
   m_TempFactor = 0.9;
   m_MeltCount = 0;
   m_TransCount = 0;
   m_EquilCount = 0;
   m_NumAborts = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumUphill = 0;
   m_NumDownhill = 0;

   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   numParams = pGroup->GetNumParams();   

   NEW_PRINT("ModelBackup", 1);
   m_pTransBackup = new ModelBackup(m_pModel);
   MEM_CHECK(m_pTransBackup);
   
   NEW_PRINT("double", numParams);
   m_pBest = new double[numParams];
   MEM_CHECK(m_pBest);

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("SAAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginSimulatedAlg", pFileName) == true)
   {
      FindToken(inFile, "EndSimulatedAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginSimulatedAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndSimulatedAlg") == NULL)
      {
         if(strstr(line, "NumInitialTrials") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumMelts);
         }
         else if(strstr(line, "TemperatureScaleFactor") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_TempFactor);
         }
         else if(strstr(line, "FinalTemperature") != NULL)
         {
            sscanf(line, "%s %s", tmp, tfnl);
            MyStrLwr(tfnl);
            if((strcmp(tfnl, "computed-vanderbilt") == 0) ||
               (strcmp(tfnl, "computed") == 0))
            {
               m_FinalTempMethod = TMETHD_VNDR;
            }
            else if(strcmp(tfnl, "computed-ben-ameur") == 0)
            {
               m_FinalTempMethod = TMETHD_BAMR;
            }
            else 
            {
               m_FinalTempMethod = TMETHD_USER;
               m_FinalTemp = atof(tfnl);
            }            
         }
         else if(strstr(line, "TransitionMethod") != NULL)
         {
            sscanf(line, "%s %s", tmp, trans);
            MyStrLwr(trans);
            if(strcmp(trans, "uniform") == 0)
            {
               m_TransitionMethod = TRANS_UNFRM;
            }
            else if(strcmp(trans, "gauss") == 0)
            {
               m_TransitionMethod = TRANS_GAUSS;
            }
         }
         else if(strstr(line, "OuterIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxOuter);
         }
         else if(strstr(line, "InnerIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxInner);
         }
         else if(strstr(line, "ConvergenceVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_StopVal);
         }
         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default algorithm setup.");
   }/* end else() */

   fclose(inFile);

   IncCtorCount();
}/* end default CTOR */

/******************************************************************************
Destroy()

Frees up memory used by the algorithm.
******************************************************************************/
void SAAlgorithm::Destroy(void)
{
   delete [] m_pBest;
   delete [] m_pMelts;
   delete m_pStats;
   delete m_pTransBackup;

   delete [] m_Finner;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void SAAlgorithm::WriteMetrics(FILE * pFile)
{
   double sd, p;
   sd = ((m_InitTemp/100.00)-m_dEavg)/3.00;
   p = exp(-m_dEavg/m_CurTemp);
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Simulated Annealing for Continuous Parameters\n");
   fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
   fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
   fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
   fprintf(pFile, "Actual Outer Iterations : %d\n", m_NumOuter);
   fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
   fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor*100.0);
   fprintf(pFile, "Initial Temperature     : %E\n", m_InitTemp);
   fprintf(pFile, "Avg. Energy Change      : %E\n", m_dEavg);
   fprintf(pFile, "Std. Dev. Energy Change : %E\n", sd);
   fprintf(pFile, "Final Temperature       : %E\n", m_CurTemp);
   fprintf(pFile, "Initial Pr[Acc]         : %.2lf%%\n", m_InitProb*100.0);
   fprintf(pFile, "Actual Final Pr[Acc]    : %.2lf%%\n", m_CurProb*100.0);
   fprintf(pFile, "Expected Final Pr[Acc]  : %.2lf%%\n", p*100.0);
   fprintf(pFile, "Melting Evals           : %d\n", m_MeltCount);
   fprintf(pFile, "Transition Evals        : %d\n", m_TransCount);
   fprintf(pFile, "Equilibration Evals     : %d\n", m_EquilCount);
   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());      
   fprintf(pFile, "Rejected Transitions    : %d\n", m_NumAborts);
   fprintf(pFile, "Uphill Transitions      : %d\n", m_NumUphill);
   fprintf(pFile, "Downhill Transitions    : %d\n", m_NumDownhill);
   fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

   m_pModel->WriteMetrics(pFile);
   if(m_CurStop <= m_StopVal)
   {
      fprintf(pFile, "Algorithm successfully converged on a solution\n");
   }
   else
   {
      fprintf(pFile, "Algorithm failed to converge on a solution, more outer iterations may be needed\n");
   }      
}/* end WriteMetrics() */

/******************************************************************************
StoreBest()

Saves the currently active parameter set into the m_pBest array.
******************************************************************************/
void SAAlgorithm::StoreBest(void)
{
   m_pModel->GetParamGroupPtr()->ReadParams(m_pBest);
} /* end StoreBest() */

/******************************************************************************
RestoreBest()

Copies the parameter set stored in the m_pBest array into the model Parameter
Group, then the model is rerun, so that all constraints, response vars and
observations are consistent.
******************************************************************************/
void SAAlgorithm::RestoreBest(void)
{
   m_pModel->GetParamGroupPtr()->WriteParams(m_pBest);
   m_pModel->Execute();
}/* end RestoreBest() */

/******************************************************************************
Optimize()

Optimize the objective function using the SA algorithm.
******************************************************************************/
void SAAlgorithm::Optimize(void)
{
   int rank, nprocs;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if(nprocs == 1)
   {
      OptimizeSerial();
   }
   else
   {
      OptimizeParallel(rank, nprocs);
   }
}/* end Optimize() */

/******************************************************************************
OptimizeSerial()

Optimize the objective function using the SA algorithm.
******************************************************************************/
void SAAlgorithm::OptimizeSerial(void)
{   
   StatusStruct pStatus;
   double curVal;
   int i;

   //write setup
   WriteSetup(m_pModel, "Simulated Annealing for Continuous Parameters");

   //read in best result from previous run, if desired
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }

   m_NumOuter = 0;

   curVal = m_pModel->Execute();
   StoreBest();
   m_MeltCount++;

   curVal = Melt(curVal);

   //write banner and initial result
   WriteBanner(m_pModel, "iter  obj. function  ", "Convergence Value");
   WriteRecord(m_pModel, 0, curVal, m_CurStop);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxOuter;
   pStatus.pct = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   for(i = 0; i < m_MaxOuter; i++)
   {
      if(IsQuit() == true){ break;}

      curVal = Equilibrate(curVal);
      m_CurTemp *= m_TempFactor;

      //write iteration result
      WriteRecord(m_pModel, (i+1), curVal, m_CurStop);
      pStatus.curIter = m_NumOuter = i+1;
      pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_MaxOuter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //check for convergence
      if(m_CurStop <= m_StopVal){pStatus.pct = 100.00; break;}

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   } /* end while() */

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write optimal results 
   WriteOptimal(m_pModel, curVal);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   //write algorithm metrics
   WriteAlgMetrics(this);
}/* end OptimizeSerial() */

/******************************************************************************
OptimizeParallel()

Optimize the objective function using the SA algorithm and multiple cpus.
******************************************************************************/
void SAAlgorithm::OptimizeParallel(int rank, int nprocs)
{
   int i;
   StatusStruct pStatus;
   double fbest;

   m_NumOuter = 0;

   if(rank == 0)
   {
      //write setup
      WriteSetup(m_pModel, "Simulated Annealing for Continuous Parameters");

      fbest = InitMaster(nprocs);

      fbest = MeltMaster(fbest, nprocs);

      //write banner and initial result
      WriteBanner(m_pModel, "iter  obj. function  ", "Convergence Value");
      WriteRecord(m_pModel, 0, fbest, m_CurStop);
      pStatus.curIter = 0;
      pStatus.maxIter = m_MaxOuter;
      pStatus.pct = 0.00;
      pStatus.numRuns = m_NumMelts + nprocs;
      WriteStatus(&pStatus);

      for(i = 0; i < m_MaxOuter; i++)
      {
         if(IsQuit() == true){ break;}

         fbest = EquilibrateMaster(fbest, nprocs);
         m_CurTemp *= m_TempFactor;

         //write iteration result
         WriteRecord(m_pModel, (i+1), fbest, m_CurStop);
         pStatus.curIter = m_NumOuter = i+1;
         if(m_CurStop <= m_StopVal)
         {
            pStatus.pct = 100.00; 
         }
         else
         {
            pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_MaxOuter;
         }
         pStatus.numRuns += m_MaxInner;
         WriteStatus(&pStatus);

         //send status update to slaves
         MPI_Bcast(&(pStatus.pct), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         //check for convergence
         if(pStatus.pct >= 100.00)
         {
            break;
         }

         //perform intermediate bookkeeping
         m_pModel->Bookkeep(false);
      } /* end for() */
   }/* end if(rank == 0) */
   else
   {
      InitSlave(rank, nprocs);
      MeltSlave(rank, nprocs);
      for(i = 0; i < m_MaxOuter; i++)
      {
         m_NumOuter = i + 1;
         if(IsQuit() == true){ break;}
         EquilibrateSlave(rank, nprocs);
         //get status update from master
         MPI_Bcast(&(pStatus.pct), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         //check for convergence
         if(pStatus.pct >= 100.00)
         {
            break;
         }
      }/* end for() */
   }/* end else() */

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(rank == 0)
   {
      //write optimal results 
      WriteOptimal(m_pModel, fbest);
      WriteStatus(&pStatus);
      //write algorithm metrics
      WriteAlgMetrics(this);
   }
}/* end OptimizeParallel() */

/******************************************************************************
InitMaster()

Master portion of initialization.
******************************************************************************/
double SAAlgorithm::InitMaster(int nprocs)
{
   MPI_Status mpi_status;
   bool bDone = false;
   int i, j, np, num_recv;
   int nslaves, sid, nxtsid;
   bool bSynch = SynchReceives();
   double f, range, lwr, upr, r, fbest;
   double * fplus;
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();
   np = pGroup->GetNumParams();

   //read in best result from previous run, if desired
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
   }

   //save initial parameters
   pGroup->ReadParams(m_pBest);

   fplus = new double[np+1];

   //output banner
   WriteInnerEval(WRITE_SMP, nprocs-1, '.');

   //assign initial work to slaves
   for(i = 1; i < nprocs; i++)
   {
      //perturb initial parameters by +/-5%
      for(j = 0; j < np; j++)
      {
         //bounds
         lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
         upr = pGroup->GetParamPtr(j)->GetUprBnd();
         //10% of range
         range = 0.1*(upr - lwr);
         //random number from -0.5 to +0.5
         r = (UniformRandom() - 0.5);
         //-5% to +5% of range
         r *= range;
         //tack on initial guess
         r += pGroup->GetParamPtr(j)->GetEstVal();
         //enforce bounds
         if(r > upr) r = upr;
         if(r < lwr) r = lwr;

         fplus[j] =  r;
      }/* end for() */
      pGroup->WriteParams(fplus);
      m_pModel->PerformParameterCorrections();

      // send work to slave
      MPI_Send(&(fplus[0]), np, MPI_DOUBLE, i, MPI_DATA_TAG, MPI_COMM_WORLD);
   }/* end for() */

   //master runs the user-specified initial value
   pGroup->WriteParams(m_pBest);
   m_pModel->PerformParameterCorrections();
   fbest = m_pModel->Execute();

   num_recv = 0;
   sid = 0;
   nxtsid = 0;
   nslaves = nprocs - 1;
   while(bDone == false)
   {
      if(bSynch == true)
      {
         sid = nxtsid + 1;
         nxtsid = (nxtsid + 1) % nslaves;
      }
      else
      {
         sid = MPI_ANY_SOURCE;
      }

      //receive result from slave and process
      MPI_Recv(&(fplus[0]), np+1, MPI_DOUBLE, sid, MPI_RESULTS_TAG, MPI_COMM_WORLD, &mpi_status);
      f = fplus[np];
      if(f < fbest)
      {
         fbest = f;
         pGroup->WriteParams(fplus);
         pGroup->ReadParams(m_pBest);
         m_pModel->SaveBest(mpi_status.MPI_SOURCE);
      }
      num_recv++;
      WriteInnerEval(num_recv, nprocs, '.');

      if(num_recv == (nprocs-1)) //quit when all slaves have reported
      {
         WriteInnerEval(WRITE_ENDED, nprocs, '.');
         bDone = true;
      }
   }/* end while() */

   //synch up processors
   MPI_Barrier(MPI_COMM_WORLD);

   delete [] fplus;

   return fbest;
} /* end InitMaster() */

/******************************************************************************
InitSlave()

Slave portion of initialization.
******************************************************************************/
void SAAlgorithm::InitSlave(int rank, int nprocs)
{
   MPI_Status mpi_status;
   int np;
   double * fplus;
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();
   np = pGroup->GetNumParams();
   
   fplus = new double[np+1];

   // receive work
   MPI_Recv(&(fplus[0]), np, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
   
   //run the model
   pGroup->WriteParams(fplus);
   fplus[np] = m_pModel->Execute();

   //send results to master
   MPI_Send(&(fplus[0]), np+1, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);

   //synch up processors
   MPI_Barrier(MPI_COMM_WORLD);

   delete [] fplus;
} /* end InitSlave() */

/******************************************************************************
Melt()

'Melts' the design space to determine the initial temperature. The initial 
temperature is computed so that a statistically large energy increase (dE) from a 
sample of random moves (melting trials) will be accepted with ~100% probability. 
That is:
   T0 = -dEmax/ln(0.99) ~ 100*dEmax
Where dEmax is computed as 2 std. deviations from the mean dE estimated from the
melting trials.

Returns objective function value at starting location.
******************************************************************************/
double SAAlgorithm::Melt(double initVal)
{
   int i; 
   double Ebest, Eprev, Ecur; //objective function values
   double Emed;               //melting trial metrics
   double dE, dEmax, dEavg;   //energy change (delta E)
   double * pdE;              //array of positive energy changes
   char c;

   //allocate space for melts (needed for convergence value)
   if(m_pMelts == NULL)
   {
      NEW_PRINT("double", m_NumMelts); 
      m_pMelts = new double[m_NumMelts];
      pdE      = new double[m_NumMelts];
      MEM_CHECK(pdE);
   }

   while(CheckOverflow(initVal*initVal) == true)
   {
     //make a random move
     GenerateRandomMove();
	  initVal = m_pModel->Execute();
   }
   
   Ebest = Ecur = initVal;
   dE = dEavg = 0.00;
   
   WriteMelt(0, m_NumMelts, '.');
   for(i = 0; i < m_NumMelts; i++)
   {
      Eprev = Ecur;
      do
      {
        //make a random move
        GenerateRandomMove();
	     Ecur  = m_pModel->Execute();
      }while(CheckOverflow(Ecur*Ecur) == true);

      dE = Ecur - Eprev;      
      m_MeltCount++;

      //accumulate values
      pdE[i] = fabs(dE);
      dEavg += pdE[i];
      m_pMelts[i] = Ecur;

      c = '+';      
      if(dE < 0.00) //decrease in energy
      {         
         c = '-';
         if(Ecur < Ebest) //new best solution?
         {
            StoreBest();
            Ebest = Ecur;
         }
      } /* end if() */      

      WriteMelt(i+1, m_MaxInner, c);
   } /* end while() */
   WriteMelt(-1, -1, '.');

   /* ------------------------------------------------------
   average or median energy change, use whichever is smaller
   to avoid influence of outliers
   ------------------------------------------------------- */
   double dEmed;
   dEavg /= (double)m_NumMelts;
   dEmed = CalcMedian(pdE, m_NumMelts);
   if(dEmed < dEavg) dEavg = dEmed;

   //compute current convergence value
   Emed = CalcMedian(m_pMelts, m_NumMelts);   
   m_CurStop = fabs((Emed - Ebest)/Emed);

   RestoreBest();   
   m_MeltCount++;

   //assign initial temperatue
   dEmax = dEavg + 3.00*CalcStdDev(pdE, m_NumMelts,CENTRAL_TEND_PCTILE);
   m_CurTemp = m_InitTemp = 100.00*dEmax;
   m_dEavg = dEavg; //store average energy change as a metric   

   //if user supplied final temp, compute corresponding reduction factor
   if(m_FinalTempMethod == TMETHD_USER)
   {
       m_TempFactor = pow((m_FinalTemp/m_InitTemp), 1.00/(double)m_MaxOuter);
   }
   //if user wants pre-computed initial and final temperatures
   else if(m_FinalTempMethod == TMETHD_BAMR)
   {
      /* -----------------------------------------------------------------
      Use algorithm described by Walid Ben-Ameur in "Computing the Initial
      Temperature of Simulated Annealing", 2004, Computational Optimization
      and Applications, vol. 29, pg. 369-385.
      ----------------------------------------------------------------- */
      bool bFirstPass = true;
      //Step 1a -- number of transitions is one less than melting trials
      int S = (m_NumMelts - 1);
      //Step 1b -- list of positive transitions from melting trials
      double * pS = m_pMelts;
      //Step 1c -- Set T1 equal to previously computed init. temperature
      double T1 = m_InitTemp;
      //Step 2 -- iteratively compute Pn (the acceptance prob.) and Tn
      double numer, denom, Pn, P0, Tn, Emax, Emin, * pT, p, dpLast, dP;
      P0 = 0.99; //initial acceptance probability
      pT = &m_InitTemp;
step2:
      dpLast = 2.00;
      p = 1.00;
      Tn = T1;
      for(int t = 0; t < 10000; t++)
      {
         numer = denom = 0.00;
         for (i = 0; i < S; i++)
         {
            if(pS[i] > pS[i+1]){ Emax = pS[i]; Emin = pS[i+1];}
            else               { Emin = pS[i]; Emax = pS[i+1];}        
            numer += exp(-Emax/Tn);
            denom += exp(-Emin/Tn);
         }
         Pn = numer/denom;
         dP = fabs(Pn - P0);
         if (dP <= 0.001) break;
         //divergence check
         if(dP >= dpLast)
         {
           p *= 2.00;
         }
         else
         {
           Tn = Tn*pow(log(Pn)/log(P0), (1.00/p));
           dpLast = dP;
         }
      }/* end for() */
      *pT = Tn;
      P0 = 0.01; //final acceptance probability
      pT = &m_FinalTemp;
      if(bFirstPass == true)
      {
         bFirstPass = false;
         goto step2; //repeat step 2 for final temperature
      }
      m_TempFactor = pow((m_FinalTemp/m_InitTemp), 1.00/(double)m_MaxOuter);
      m_CurTemp = m_InitTemp;
   }/* end else if() */
   else if(m_FinalTempMethod == TMETHD_VNDR)
   {
      m_InitTemp = -m_dEavg/log(0.99); 
      m_FinalTemp = -m_dEavg/log(0.01);
      m_TempFactor = pow((m_FinalTemp/m_InitTemp), 1.00/(double)m_MaxOuter);
      m_CurTemp = m_InitTemp;
   }/* end else if() */

   //if reduction rate is not valid, compute internally
   if((m_TempFactor >= 1.00) || (m_TempFactor <= 0.00))
   {
      LogError(ERR_BAD_ARGS, "Invalid temperature reduction rate; using internally calculated value");
      m_TempFactor = pow((1.00/m_InitTemp), 1.00/(double)m_MaxOuter);
   }
   
   delete [] pdE;

   return (m_pModel->GetObjFuncVal());
}/* end Melt() */

/******************************************************************************
MeltMaster()

Supervise the melting operation.
******************************************************************************/
double SAAlgorithm::MeltMaster(double fbest, int nprocs)
{
   MPI_Status mpi_status;
   int i, signal, np, nstops, num_recv, sid, nxtsid, nslaves; 
   double * fplus;
   double * pdE;              //array of positive energy changes
   double fcur, fprev;
   double Emed, dE, dEmax, dEavg;   //energy change (delta E)
   char c;
   bool bDone = false;
   bool bSynch = SynchReceives();
   ParameterGroup * pGroup;

   //allocate space for melts
   if(m_pMelts == NULL)
   {
      NEW_PRINT("double", m_NumMelts); 
      m_pMelts = new double[m_NumMelts];
      pdE = new double[m_NumMelts];
      MEM_CHECK(pdE);
   }

   pGroup = m_pModel->GetParamGroupPtr();
   np = pGroup->GetNumParams();
   
   fplus = new double[np+1];
   fcur = fprev = fbest;

   WriteMelt(0, m_NumMelts, '.');
   nstops = 0;

   //send initial set of melts to slaves
   for(i = 1; i < nprocs; i++)
   {
      if(i <= m_NumMelts)
      {
         //make a random move
         GenerateRandomMove();
         m_pModel->PerformParameterCorrections();
         pGroup->ReadParams(fplus);

         // send move to slave
         signal = APSA_DO_WORK;
         MPI_Send(&signal,1,MPI_INT,i,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
         MPI_Send(&(fplus[0]), np, MPI_DOUBLE, i, MPI_DATA_TAG, MPI_COMM_WORLD);
      }
      else
      {
         signal = APSA_STOP_WORK;
         MPI_Send(&signal,1,MPI_INT,i,MPI_REQUEST_TAG,MPI_COMM_WORLD);
         nstops++;
      }
   }/* end for(each slave) */
   
   //receive results from slaves and send more work as needed
   num_recv = 0;
   dE = dEavg = 0.00;
   nxtsid = 0;
   sid = 0;
   nslaves = nprocs - 1;
   while(bDone == false)
   {

      if(bSynch == true)
      {
         sid = nxtsid + 1;
         nxtsid = (nxtsid + 1) % nslaves;
      }
      else
      {
         sid = MPI_ANY_SOURCE;
      }

      //receive result from slave and process
      MPI_Recv(&(fplus[0]), np+1, MPI_DOUBLE, sid, MPI_RESULTS_TAG, MPI_COMM_WORLD, &mpi_status);     
      sid = mpi_status.MPI_SOURCE;
      fcur = fplus[np];

      m_MeltCount++;

      //archive values
      if(CheckOverflow(fcur*fcur) == false)
      {
         dE = fcur - fprev;
         printf("fcur = %E, fprev = %E, dE = %E\n", fcur, fprev, dE);
         fprev = fcur;
       
         if(num_recv < m_NumMelts)
         {
            m_pMelts[num_recv] = fcur;
            //accumulate values
            pdE[num_recv] = fabs(dE);
            dEavg += pdE[num_recv];
         }

         //update best solution, if new one is found
         c = '+';
         if(fcur < fbest)
         {         
            c = '-';
            pGroup->WriteParams(fplus);
            pGroup->ReadParams(m_pBest);
            m_pModel->SaveBest(mpi_status.MPI_SOURCE);
            fbest = fcur;
         } /* end if() */

         WriteMelt(num_recv+1, m_MaxInner, c);

         num_recv++;
      }
      else //discard result
      {
         i--;
      }

      //assign more work
      if(i <= m_NumMelts)
      {
         //make a random move
         GenerateRandomMove();
         pGroup->ReadParams(fplus);
         m_pModel->PerformParameterCorrections();

         // send move to slave
         signal = APSA_DO_WORK;
         MPI_Send(&signal,1,MPI_INT,sid,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
         MPI_Send(&(fplus[0]), np, MPI_DOUBLE, sid, MPI_DATA_TAG, MPI_COMM_WORLD);
         i++;
      }
      else // send stop work message to the slave
      {
         signal = APSA_STOP_WORK;
         MPI_Send(&signal,1,MPI_INT,sid,MPI_REQUEST_TAG,MPI_COMM_WORLD);
         nstops++;
         if(nstops == (nprocs-1)) //quit when all slaves have been told to stop
         {
            //WriteInnerEval(WRITE_ENDED, m_NumMelts, '.');
            bDone = true;
         }
      }
   }/* end while() */

   //done melting
   WriteMelt(-1, -1, '.');

/*
   for(i = 0; i < m_NumMelts; i++)
   {
     printf("Melt[%d] = %E\n",i, m_pMelts[i]);
     fflush(stdout);
   }
*/

   /* ------------------------------------------------------
   average or median energy change, use whichever is smaller
   to avoid influence of outliers
   ------------------------------------------------------- */
   double dEmed;
   dEavg /= (double)m_NumMelts;
   dEmed = CalcMedian(pdE, m_NumMelts);
   if(dEmed < dEavg) dEavg = dEmed;

   //compute current convergence value
   Emed = CalcMedian(m_pMelts, m_NumMelts);   
   m_CurStop = fabs((Emed - fbest)/Emed);

   pGroup->WriteParams(m_pBest);

   //assign initial temperatue --- shoot for 95th percentile of energy change
   dEmax = dEavg + 2.00*CalcStdDev(pdE, m_NumMelts,CENTRAL_TEND_PCTILE);
   m_CurTemp = m_InitTemp = 100.00*dEmax;
   m_dEavg = dEavg; //store average energy change as a metric   

   //if user supplied final temp, compute corresponding reduction factor
   if(m_FinalTempMethod == TMETHD_USER)
   {
       m_TempFactor = pow((m_FinalTemp/m_InitTemp), 1.00/(double)m_MaxOuter);
   }
   //if user wants pre-computed initial and final temperatures
   else if(m_FinalTempMethod == TMETHD_BAMR)
   {
      /* -----------------------------------------------------------------
      Use algorithm described by Walid Ben-Ameur in "Computing the Initial
      Temperature of Simulated Annealing", 2004, Computational Optimization
      and Applications, vol. 29, pg. 369-385.
      ----------------------------------------------------------------- */
      bool bFirstPass = true;
      //Step 1a -- number of transitions is one less than melting trials
      int S = (m_NumMelts - 1);
      //Step 1b -- list of positive transitions from melting trials
      double * pS = m_pMelts;
      //Step 1c -- Set T1 equal to previously computed init. temperature
      double T1 = m_InitTemp;
      //Step 2 -- iteratively compute Pn (the acceptance prob.) and Tn
      double numer, denom, Pn, P0, Tn, Emax, Emin, * pT, p, dpLast, dP;
      P0 = 0.99; //initial acceptance probability
      pT = &m_InitTemp;
step2:
      dpLast = 2.00;
      p = 1.00;
      Tn = T1;
      for(int t = 0; t < 10000; t++)
      {
         numer = denom = 0.00;
         for (i = 0; i < S; i++)
         {
            if(pS[i] > pS[i+1]){ Emax = pS[i]; Emin = pS[i+1];}
            else               { Emin = pS[i]; Emax = pS[i+1];}        
            numer += exp(-Emax/Tn);
            denom += exp(-Emin/Tn);
         }
         Pn = numer/denom;
         dP = fabs(Pn - P0);
         if (dP <= 0.001) break;
         //divergence check
         if(dP >= dpLast)
         {
           p *= 2.00;
         }
         else
         {
           Tn = Tn*pow(log(Pn)/log(P0), (1.00/p));
           dpLast = dP;
         }
      }/* end for() */
      *pT = Tn;
      P0 = 0.01; //final acceptance probability
      pT = &m_FinalTemp;
      if(bFirstPass == true)
      {
         bFirstPass = false;
         goto step2; //repeat step 2 for final temperature
      }
      m_TempFactor = pow((m_FinalTemp/m_InitTemp), 1.00/(double)m_MaxOuter);
      m_CurTemp = m_InitTemp;
   }/* end else if() */
   else if(m_FinalTempMethod == TMETHD_VNDR)
   {
      m_InitTemp = -m_dEavg/log(0.99); 
      m_FinalTemp = -m_dEavg/log(0.01);
      m_TempFactor = pow((m_FinalTemp/m_InitTemp), 1.00/(double)m_MaxOuter);
      m_CurTemp = m_InitTemp;
   }/* end else if() */

   //if reduction rate is not valid, compute internally
   if((m_TempFactor >= 1.00) || (m_TempFactor <= 0.00))
   {
      LogError(ERR_BAD_ARGS, "Invalid temperature reduction rate; using internally calculated value");
      m_TempFactor = pow((1.00/m_InitTemp), 1.00/(double)m_MaxOuter);
   }
   
   //synch up processors
   MPI_Barrier(MPI_COMM_WORLD);

   delete [] pdE;
   delete [] fplus;
   return fbest;
}/* end MeltMaster() */

/******************************************************************************
MeltSlave()

Perform model evaluations assigned by the master.
******************************************************************************/
void SAAlgorithm::MeltSlave(int rank, int nprocs)
{
   int np;
   double * fplus;
   ParameterGroup * pGroup;
   MPI_Status mpi_status;
   int signal;
   bool bDone = false;

   pGroup = m_pModel->GetParamGroupPtr();
   np = pGroup->GetNumParams();   
   fplus = new double[np+1];

   while(bDone == false)
   {
      MPI_Recv(&signal,1,MPI_INT,0,MPI_REQUEST_TAG,MPI_COMM_WORLD, &mpi_status); 
      if(signal == APSA_DO_WORK)
      {
         MPI_Recv(&(fplus[0]), np, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
         pGroup->WriteParams(fplus);
         fplus[np] = m_pModel->Execute();
         MPI_Send(&(fplus[0]), np+1, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);
      }/* end if() */
      else
      {
         bDone = true;
      }
   }/* end while() */

   delete [] fplus;

   //synch up processors
   MPI_Barrier(MPI_COMM_WORLD);
}/* end MeltSlave() */

/******************************************************************************
Equilibrate()

Allows the system to come to equilibrium at the current temperature. This is 
accomplished by making a number of moves equal to m_MaxInner and choosing 
the best one. Returns the objective function value of the best move.
******************************************************************************/
double SAAlgorithm::Equilibrate(double initVal)
{
   int m,n; 
   double bestVal, curVal, lastVal, avgVal, median;
   char c;
   ParameterGroup * pGroup;
 
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();

   if(m_Finner == NULL)
   {
      //storage for inner loop f(x) calculations (used in convergence test)
      NEW_PRINT("double", m_MaxInner);
      m_Finner = new double [m_MaxInner];
      MEM_CHECK(m_Finner);
   }/* end if() */

   bestVal = curVal = initVal;

   WriteInnerEval(WRITE_SA, m_MaxInner, '.');

   //initialize probability metrics for this set of tranisitions
   m_NumProbTests = 0;
   m_TotProb      = 0.00;
   avgVal = 0.00;

   for(m = 0; m < m_MaxInner; m++)
   {
      lastVal = curVal;
      if(m_TransitionMethod == TRANS_UNFRM)
        curVal = Transition(curVal);
      else 
        curVal = GaussTransition(curVal);

      avgVal += curVal;
      m_Finner[m] = curVal;

      //update current best
      if(curVal < bestVal)
      {         
         StoreBest();
         bestVal = curVal;
      } /* end if() */

      if(curVal < lastVal)      { c = '-';}
      else if(curVal == lastVal){ c = '.';}
      else                      { c = '+';}
      WriteInnerEval(m+1, m_MaxInner, c);
   } /* end for() */

   //compute current convergence value (Eqn 18 of paper)
   avgVal /= (double)m_MaxInner;
   median = CalcMedian(m_Finner, m_MaxInner);
   m_CurStop = fabs((median - bestVal)/median);

   //compute avg. probability of acceptance for the equilibration
   if(m_NumProbTests > 0)
   { 
      m_CurProb = m_TotProb / (double)m_NumProbTests;
      if(m_InitProb < 0.00){ m_InitProb = m_CurProb;}
   }

   WriteInnerEval(WRITE_ENDED, m_MaxInner, '.');

   RestoreBest();
   m_EquilCount++;

   return bestVal;
} /* end Equilibrate() */

/******************************************************************************
EquilibrateMaster()

Allows the system to come to equilibrium at the current temperature. This is 
accomplished by making a number of moves equal to m_MaxInner and choosing 
the best one. Returns the objective function value of the best move.
******************************************************************************/
double SAAlgorithm::EquilibrateMaster(double fbest, int nprocs)
{
   MPI_Status mpi_status;
   int m, n, num_recv, signal, nstops; 
   double fcur, flast,  median, finit;
   int sid, nxtsid, nslaves;
   char c;
   bool bDone = false;
   bool bSynch = SynchReceives();
   ParameterGroup * pGroup;
 
   //allocate m_x array, if needed
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   if(m_Finner == NULL)
   {
      //storage for inner loop f(x) calculations
      NEW_PRINT("double", m_MaxInner);
      m_Finner = new double [m_MaxInner];
      MEM_CHECK(m_Finner);
   }/* end if() */

   fcur = finit = fbest;

   WriteInnerEval(WRITE_SA, m_MaxInner, '.');

   //initialize probability metrics for this set of tranisitions
   m_NumProbTests = 0;
   m_TotProb = 0.00;

   //prime the pump --- get slaves started on initial set of transitions
   for(m = 1; m < nprocs; m++)
   {
      if(m_TransitionMethod == TRANS_UNFRM)
      {
         TransitionSend(fcur, m);
      }
      else
      {
         GaussTransitionSend(fcur, m);
      }
   }/* end for() */

   nstops = 0;
   num_recv = 0;
   sid = nxtsid = 0;
   nslaves = nprocs - 1;
   while(bDone == false)
   {
      flast = fcur;

      if(bSynch == true)
      {
         sid = nxtsid + 1;
         nxtsid = (nxtsid + 1) % nslaves;
      }
      else
      {
         sid = MPI_ANY_SOURCE;
      }
      mpi_status.MPI_SOURCE = sid;

      if(m_TransitionMethod == TRANS_UNFRM)
      {
         fcur = TransitionRecv(fcur, &(mpi_status));
      }
      else
      {
         fcur = GaussTransitionRecv(fcur, &(mpi_status));
      }

      m_Finner[num_recv] = fcur;

      //update current best
      if(fcur < fbest)
      {         
         pGroup->ReadParams(m_pBest);
         fbest = fcur;
         m_pModel->SaveBest(mpi_status.MPI_SOURCE);
      } /* end if() */

      if(fcur < flast) { c = '-';}
      else if(fcur == flast){ c = '.';}
      else { c = '+';}
      WriteInnerEval(num_recv+1, m_MaxInner, c);
    
      num_recv++;
  
      //assign more work
      if(m <= m_MaxInner)
      {
         if(m_TransitionMethod == TRANS_UNFRM)
         {
            TransitionSend(fcur, mpi_status.MPI_SOURCE);
         }
         else
         {
            GaussTransitionSend(fcur, mpi_status.MPI_SOURCE);
         }
         m++;
      }/* end if() */
      else // send stop work message to the slave
      {
         signal = APSA_STOP_WORK;
         MPI_Send(&signal,1,MPI_INT,mpi_status.MPI_SOURCE,MPI_REQUEST_TAG,MPI_COMM_WORLD);
         nstops++;
         if(nstops == (nprocs-1)) //quit when all slaves have been told to stop
         {
            //WriteInnerEval(WRITE_ENDED, m_MaxInner, '.');
            bDone = true;
         }
      }/* end else() */
   }/* end while() */

   //compute current convergence value   
   median = CalcMedian(m_Finner, m_MaxInner);
   m_CurStop = fabs((median - fbest)/median);

   //compute avg. probability of acceptance for the equilibration
   if(m_NumProbTests > 0)
   { 
      m_CurProb = m_TotProb / (double)m_NumProbTests;
      if(m_InitProb < 0.00){ m_InitProb = m_CurProb;}
   }

   WriteInnerEval(WRITE_ENDED, m_MaxInner, '.');

   pGroup->WriteParams(m_pBest);;
   m_EquilCount++;

   return fbest;
}/* end EquilibrateMaster() */

/******************************************************************************
EquilibrateSlave()

Peform model evaluation on behalf of master processor.
******************************************************************************/
void SAAlgorithm::EquilibrateSlave(int rank, int nprocs)
{   
   double * fplus;
   MPI_Status mpi_status;
   int signal, np;

   np = m_pModel->GetParamGroupPtr()->GetNumParams();
   fplus = new double[np+1];
   signal = APSA_DO_WORK; 
   while(signal == APSA_DO_WORK)
   {
      MPI_Recv(&signal,1,MPI_INT,0,MPI_REQUEST_TAG,MPI_COMM_WORLD, &mpi_status); 
      if(signal == APSA_DO_WORK)
      {
         MPI_Recv(&(fplus[0]), np, MPI_DOUBLE, 0, MPI_DATA_TAG, MPI_COMM_WORLD, &mpi_status);
         m_pModel->GetParamGroupPtr()->WriteParams(fplus);
	      fplus[np] = m_pModel->Execute();
         MPI_Send(&(fplus[0]), np+1, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);
      }/* end if() */
   }/* end while() */
}/* end EquilibrateSlave() */

/******************************************************************************
Transition()

Attempts to make a move from the current paramter set (location). For each 
move attempted, the resulting objective function value is tested against the
acceptance criteria (either the move reduces ob. func. or a randomly generated
number is less than the acceptance probability. Returns the value of the obj.
function at the revised location.

The amount of parameter change induced by a transition is defined by the size
of the 'neighborhood' of adjacent solutions, which for continuous
SA is limited to +/-10% of the parameter range.
******************************************************************************/
double SAAlgorithm::Transition(double initVal)
{
   int i, n;
   double increase, prob, r, fcur;
   double upr, lwr, range, val, curval;
   double range10pct, vavg, vmin, vmax;
   double a;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   a = (double)(m_pModel->GetCounter() - m_NumMelts) / (double) (m_MaxOuter*m_MaxInner);

	//store initial Parameter values
   m_pTransBackup->Store();
  
   /*-------------------------------------------------------------
   Make a move --- randomly perturb parameters relative to range
   or relative to current value 
   -------------------------------------------------------------*/	
   unsigned int rangeType = MyRand()%2;   
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   for(i = 0; i < n; i++)
   {
      //perterb the parameter to a neighboring state (+/- 10% of range)
      pParam = pGroup->GetParamPtr(i);
      curval = pParam->GetEstVal();
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      if(rangeType == 0) //absolute range
      {
         range = upr - lwr;
      }
      else //relative range
      {
         range = curval;
         if(range <= NEARLY_ZERO)
         {
            range = upr - lwr;
         }
      }
      
      if((range >= 1.00) && (range < 10.00) && 
         (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;

      /* +/-10% of the range = range/5 */
      range10pct = range/5;

      /* uniform random number between 0 and 1 */
      r = ((double)MyRand() / (double)MY_RAND_MAX);
      
      /* minimum value after perturbation --- logic handles edge cases */
      vavg = (upr + lwr)*0.5;
      if(curval < vavg)
      {
         vmin = curval - 0.5*range10pct;
         if(vmin < lwr) vmin = lwr;
         vmax = vmin + range10pct;
      }
      else
      {
         vmax = curval + 0.5*range10pct;
         if(vmax > upr) vmax = upr;
         vmin = vmax - range10pct;
      }

      /* adjusted value = (10% range * random) + minimum */
      val = (range10pct * r) + vmin;
            
      //enforce limits by moving half the distance --- this shouldn't happen given logic above!
      if(val > upr){val = (upr+curval)/2.00; m_NumUprViols++;}
      if(val < lwr){val = (curval+lwr)/2.00; m_NumLwrViols++;}
      double bst = m_pTransBackup->GetParam(i);
      val = TelescopicCorrection(lwr, upr, bst, a, val);
      pParam->SetEstVal(val);
   }
   m_pModel->PerformParameterCorrections();
	fcur = m_pModel->Execute();
   m_TransCount++;

   /*---------------------------------------
   Accept or reject move.
   ---------------------------------------*/	
	if(fcur <= initVal)
	{
      m_NumDownhill++;
   } /* end if() */
   else /* test against acceptance probability */
   {			
      increase = fcur - initVal;

      // generate a random between 0 and 1
      r = (double)MyRand() / (double)MY_RAND_MAX;

      //check if the increase is acceptable
      prob = exp(-increase / m_CurTemp);

      //probability metrics
      m_TotProb += prob;
      m_NumProbTests++;

      if(prob >= r)
      {
         // accept the move     
         m_NumUphill++;
      } /* end if() */
      else
      {
         m_pTransBackup->SemiRestore();         
         m_NumAborts++;
         fcur = initVal;
      } /* end else() */
   }/* end else() */

   return fcur;
} /* end Transition() */

/******************************************************************************
TransitionSend()

Send a transition request to a slave processor.
******************************************************************************/
void SAAlgorithm::TransitionSend(double finit, int whichProc)
{

   int i, n;
   double upr, lwr, range, range10pct, val, curval;
   double a, r, vavg, vmin, vmax;
   double * curParams, * newParams;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   a = (double)(m_TransCount) / (double) (m_MaxOuter*m_MaxInner);

   /*-------------------------------------------------------------
   Make a move --- randomly perturb parameters relative to range
   or relative to current value 
   -------------------------------------------------------------*/	
   unsigned int rangeType = MyRand()%2;   
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   curParams = new double[n];
   pGroup->ReadParams(curParams);
   newParams = new double[n];
   for(i = 0; i < n; i++)
   {
      //perterb the parameter to a neighboring state (+/- 10% of range)
      pParam = pGroup->GetParamPtr(i);
      curval = pParam->GetEstVal();
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      if(rangeType == 0) //absolute range
      {
         range = upr - lwr;
      }
      else //relative range
      {
         range = curval;
         if(range <= NEARLY_ZERO)
         {
            range = upr - lwr;
         }
      }

      if((range >= 1.00) && (range < 10.00) && 
         (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;
     
      /* +/-10% of the range = range/5 */
      range10pct = range/5;

      /* uniform random number between 0 and 1 */
      r = ((double)MyRand() / (double)MY_RAND_MAX);
      
      /* minimum value after perturbation --- logic handles edge cases */
      vavg = (upr + lwr)*0.5;
      if(curval < vavg)
      {
         vmin = curval - 0.5*range10pct;
         if(vmin < lwr) vmin = lwr;
         vmax = vmin + range10pct;
      }
      else
      {
         vmax = curval + 0.5*range10pct;
         if(vmax > upr) vmax = upr;
         vmin = vmax - range10pct;
      }

      /* adjusted value = (10% range * random) + minimum */
      val = (range10pct * r) + vmin;
            
      //enforce limits by moving half the distance --- this shouldn't happen given logic above!
      if(val > upr){printf("val = %E is too big\n", val); val = (upr+curval)/2.00; m_NumUprViols++;}
      if(val < lwr){printf("val = %E is too small\n", val); val = (curval+lwr)/2.00; m_NumLwrViols++;}
      double bst = m_pBest[i];
      val = TelescopicCorrection(lwr, upr, bst, a, val);
      pParam->SetEstVal(val);
   }
   m_pModel->PerformParameterCorrections();
   pGroup->ReadParams(newParams);

   // send work to slave
   int signal = APSA_DO_WORK;
   MPI_Send(&signal,1,MPI_INT,whichProc,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
   MPI_Send(&(newParams[0]), n, MPI_DOUBLE, whichProc, MPI_DATA_TAG, MPI_COMM_WORLD);

   m_TransCount++;

   //restore parameter values
   pGroup->WriteParams(curParams);

   delete [] curParams;
   delete [] newParams;
}/* end TransitionSend() */

/******************************************************************************
TransitionRecv()

Receive a transition result from a slave processor and process.
******************************************************************************/
double SAAlgorithm::TransitionRecv(double finit, MPI_Status * status)
{
   double * fplus;
   int n, sid;
   double increase, prob, r, fcur;
   ParameterGroup * pGroup;
  
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   fplus = new double[n+1];

   //receive result from slave and process
   sid = status->MPI_SOURCE;
   MPI_Recv(&(fplus[0]), n+1, MPI_DOUBLE, sid, MPI_RESULTS_TAG, MPI_COMM_WORLD, status);

	fcur = fplus[n];

   /*---------------------------------------
   Accept or reject move.
   ---------------------------------------*/	
	if(fcur <= finit)
	{
      m_NumDownhill++;

      //accept move
      pGroup->WriteParams(fplus);
   } /* end if() */
   else /* test against acceptance probability */
   {			
      increase = fcur - finit;

      // generate a random between 0 and 1
      r = (double)MyRand() / (double)MY_RAND_MAX;

      //check if the increase is acceptable
      prob = exp(-increase / m_CurTemp);

      //probability metrics
      m_TotProb += prob;
      m_NumProbTests++;

      if(prob >= r)
      {
         // accept the move     
         m_NumUphill++;

         pGroup->WriteParams(fplus);         
      } /* end if() */
      else
      {         
         m_NumAborts++;
         fcur = finit;
      } /* end else() */
   }/* end else() */

   delete [] fplus;

   return fcur;
}/* end TransitionRecv() */

/******************************************************************************
GaussTransition()

Attempts to make a move from the current paramter set (location). For each 
move attempted, the resulting objective function value is tested against the
acceptance criteria (either the move reduces ob. func. or a randomly generated
number is less than the acceptance probability. Returns the value of the obj.
function at the revised location.

The amount of parameter change induced by a transition is defined by the size
of the 'neighborhood' of adjacent solutions, which for continuous
SA is limited to +/-10% of the parameter range - centered on current parameter
location.
******************************************************************************/
double SAAlgorithm::GaussTransition(double initVal)
{
   int i, n;
   double increase, prob, r, a;
   double upr, lwr, val, curVal;
   double range, sdmax, sdi;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

   a = (double)(m_pModel->GetCounter() - m_NumMelts) / (double) (m_MaxOuter*m_MaxInner);

	//store initial Parameter values
   m_pTransBackup->Store();
  
   /*---------------------------------------
   Make a move.
   ---------------------------------------*/	
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   double sd = sqrt(MyMax(NEARLY_ZERO,fabs(initVal))/(double)n);

   for(i = 0; i < n; i++)
   {
      //perterb the parameter to a neighboring state (+/- 10% of range)
      pParam = pGroup->GetParamPtr(i);
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      val = curVal = pParam->GetEstVal();

      //1 std. dev. should not cover more than 68% of the range
      range = upr - lwr;
      sdmax = range * 0.68;
      if(sd > sdmax) sdi = sdmax;
      else sdi = sd;

      val = curVal = pParam->GetEstVal();

      // epsilon perturbation using normal distribution
      // centered on childVal with standard deviation estimated 
      // by the fitness value
      val = MyGaussRand(curVal, sdi);

      //enforce parameter limits
      r = 2*((double)MyRand() / (double)MY_RAND_MAX)-1.00; //-1 to 1
      if((val > upr) && (r >= 0.00)) val = curVal + (upr-curVal)*r;
      if((val > upr) && (r <  0.00)) val = curVal + (curVal-lwr)*r;
      if((val < lwr) && (r >= 0.00)) val = curVal - (curVal-lwr)*r;
      if((val < lwr) && (r <  0.00)) val = curVal - (upr-curVal)*r;

      double bst = m_pBest[i];
      val = TelescopicCorrection(lwr, upr, bst, a, val);
      pParam->SetEstVal(val);
   }
   m_pModel->PerformParameterCorrections();
	curVal = m_pModel->Execute();
   m_TransCount++;

   /*---------------------------------------
   Accept or reject move.
   ---------------------------------------*/	
	if(curVal <= initVal)
	{
      m_NumDownhill++;
   } /* end if() */
   else /* test against acceptance probability */
   {			
      increase = curVal - initVal;

      // generate a random between 0 and 1
      r = (double)MyRand() / (double)MY_RAND_MAX;

      //check if the increase is acceptable
      prob = exp(-increase / m_CurTemp);

      //probability metrics
      m_TotProb += prob;
      m_NumProbTests++;

      if(prob >= r)
      {
         // accept the move     
         m_NumUphill++;
      } /* end if() */
      else
      {
         m_pTransBackup->SemiRestore();         
         m_NumAborts++;
         curVal = initVal;
      } /* end else() */
   }/* end else() */

   return curVal;
} /* end GaussTransition() */

/******************************************************************************
GaussTransitionSend()

Send a gauss transition request to a slave processor.
******************************************************************************/
void SAAlgorithm::GaussTransitionSend(double finit, int whichProc)
{
   int i, n;
   double r, a, sdmax, sdi;
   double upr, lwr, range, val, curVal;
   double * curParams, * newParams;
   ParameterGroup * pGroup;
   ParameterABC * pParam;
  
   a = (double)(m_TransCount) / (double) (m_MaxOuter*m_MaxInner);

   /*---------------------------------------
   Make a move.
   ---------------------------------------*/	
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   curParams = new double[n];
   pGroup->ReadParams(curParams);
   newParams = new double[n];
   double sd = sqrt(MyMax(NEARLY_ZERO,fabs(finit))/(double)n);

   for(i = 0; i < n; i++)
   {
      //perterb the parameter to a neighboring state (+/- 10% of range)
      pParam = pGroup->GetParamPtr(i);
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();

      //1 std. dev. should not cover more than 68% of the range
      range = upr - lwr;
      sdmax = range * 0.68;
      if(sd > sdmax) sdi = sdmax;
      else sdi = sd;

      val = curVal = pParam->GetEstVal();

      // epsilon perturbation using normal distribution
      // centered on childVal with standard deviation estimated 
      // by the fitness value
      val = MyGaussRand(curVal, sdi);

      //enforce parameter limits
      r = 2*((double)MyRand() / (double)MY_RAND_MAX)-1.00; //-1 to 1
      if((val > upr) && (r >= 0.00)) val = curVal + (upr-curVal)*r;
      if((val > upr) && (r <  0.00)) val = curVal + (curVal-lwr)*r;
      if((val < lwr) && (r >= 0.00)) val = curVal - (curVal-lwr)*r;
      if((val < lwr) && (r <  0.00)) val = curVal - (upr-curVal)*r;

      double bst = m_pBest[i];
      val = TelescopicCorrection(lwr, upr, bst, a, val);

      pParam->SetEstVal(val);
   }/* end for() */
   m_pModel->PerformParameterCorrections();
   pGroup->ReadParams(newParams);

   // send work to slave
   int signal = APSA_DO_WORK;
   MPI_Send(&signal,1,MPI_INT,whichProc,MPI_REQUEST_TAG,MPI_COMM_WORLD); 
   MPI_Send(&(newParams[0]), n, MPI_DOUBLE, whichProc, MPI_DATA_TAG, MPI_COMM_WORLD);
	
   m_TransCount++;

   //restore parameter values
   pGroup->WriteParams(curParams);
}/* end GaussTransitionSend() */

/******************************************************************************
GaussTransitionRecv()

Receive a Gauss transition result from a slave processor and process.
******************************************************************************/
double SAAlgorithm::GaussTransitionRecv(double finit, MPI_Status * status)
{
   double * fplus;
   int n, sid;
   double increase, prob, r;
   double fcur;
   ParameterGroup * pGroup;

   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   fplus = new double[n+1];

   //receive result from slave and process
   sid = status->MPI_SOURCE;
   MPI_Recv(&(fplus[0]), n+1, MPI_DOUBLE, sid, MPI_RESULTS_TAG, MPI_COMM_WORLD, status);

	fcur = fplus[n];

   /*---------------------------------------
   Accept or reject move.
   ---------------------------------------*/	
	if(fcur <= finit)
	{
      m_NumDownhill++;

      //accept move
      pGroup->WriteParams(fplus);
   } /* end if() */
   else /* test against acceptance probability */
   {			
      increase = fcur - finit;

      // generate a random between 0 and 1
      r = (double)MyRand() / (double)MY_RAND_MAX;

      //check if the increase is acceptable
      prob = exp(-increase / m_CurTemp);

      //probability metrics
      m_TotProb += prob;
      m_NumProbTests++;

      if(prob >= r)
      {
         // accept the move     
         m_NumUphill++;

         pGroup->WriteParams(fplus);
      } /* end if() */
      else
      {
         m_NumAborts++;
         fcur = finit;
      } /* end else() */
   }/* end else() */

   delete [] fplus;

   return fcur;
}/* end GaussTransitionRecv() */

/******************************************************************************
Calibrate()

Calibrate the model using the SA algorithm.
******************************************************************************/
void SAAlgorithm::Calibrate(void)
{ 
   FILE * pFile;
   char fileName[DEF_STR_SZ];
   int id;

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

      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      m_pStats->WriteStats(stdout);
   }
} /* end Calibrate() */

/******************************************************************************
GenerateRandomMove()

Adjusts a single parameter by a random displacement.
******************************************************************************/
void SAAlgorithm::GenerateRandomMove(ParameterABC * pParam)
{   
   double r; //random number
   double range, range10pct; //width of move set
   double curVal; //current parameter value
   double upr; //upper limit of current move
   double lwr; //lower limit of current move
   double adjVal; //adjusted value 
   double vmin, vmax, vavg;

   curVal = pParam->GetEstVal();
   upr = pParam->GetUprBnd();
   lwr = pParam->GetLwrBnd();
   range = upr - lwr;
      
   if((range >= 1.00) && (range < 10.00) && 
      (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;

   /* +/-10% of the range = range/5 */
   range10pct = range/5;

   /* uniform random number between 0 and 1 */
   r = ((double)MyRand() / (double)MY_RAND_MAX);
      
   /* minimum value after perturbation --- logic handles edge cases */
   vavg = (upr + lwr)*0.5;
   if(curVal < vavg)
   {
      vmin = curVal - 0.5*range10pct;
      if(vmin < lwr) vmin = lwr;
      vmax = vmin + range10pct;
   }
   else
   {
      vmax = curVal + 0.5*range10pct;
      if(vmax > upr) vmax = upr;
      vmin = vmax - range10pct;
   }

   /* adjusted value = (10% range * random) + minimum */
   adjVal = (range10pct * r) + vmin;
            
   //enforce limits by moving half the distance --- this shouldn't happen given logic above!
   if(adjVal > upr){adjVal = (upr+curVal)/2.00; m_NumUprViols++;}
   if(adjVal < lwr){adjVal = (curVal+lwr)/2.00; m_NumLwrViols++;}
   
   pParam->SetEstVal(adjVal);
} /* end GenerateRandomMove() */

/******************************************************************************
GenerateRandomMove()

Adjusts parameters by a random displacement.
******************************************************************************/
void SAAlgorithm::GenerateRandomMove(void)
{
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   int numParams;
   int i;

   pGroup = m_pModel->GetParamGroupPtr();   
   numParams = pGroup->GetNumParams();
   for(i = 0; i < numParams; i++)
   {
      pParam = pGroup->GetParamPtr(i);
      GenerateRandomMove(pParam);
   }
} /* end GenerateRandomMove() */

/******************************************************************************
SA_Program()

Calibrate the model using the SA algorithm.
******************************************************************************/
void SA_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("SAAlgorithm", 1);
   SAAlgorithm * SA = new SAAlgorithm(model);
   MEM_CHECK(SA);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ SA->Calibrate(); }
   else { SA->Optimize(); }

   delete SA;
   delete model;
} /* end SA_Program() */


