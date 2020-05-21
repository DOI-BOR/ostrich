/******************************************************************************
File     : RejectionSampler.cpp
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Rejection Sampling Algorithm

Version History
06-23-10    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "RejectionSampler.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObservationGroup.h"
#include "SuperMUSE.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "Utility.h"
#include "SuperMuseUtility.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.

If bMCMC is true, then use Metropolis MCMC sampler, otherwise use Rejection
Sampler.
******************************************************************************/
RejectionSampler::RejectionSampler(ModelABC * pModel, bool bMCMC)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pAccepted = NULL;
   m_pSamples = NULL;
   m_MaxSamples = 0;
   m_NumDesired = 0;
   m_NumBurnIn = 0;
   m_NumFound = 0;
   m_SamplesPerIter = 0;
   m_MinWSSE = HUGE_VAL;
   m_LastWSSE = HUGE_VAL;
   m_bStedinger = true;
   m_bBeven = false;
   m_ShapeFactor = 0.5; //RMSE
   m_bMetropolis = bMCMC; //if true use Metropolis MCMC sampler

   //fraction by which to constrict parameter bounds after each iteration
   m_TelescopeRate = 1.00; 

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   m_pBuf    = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the rejection sampler and it's member variables.
******************************************************************************/
void RejectionSampler::Destroy(void)
{
   int i;
   for(i = 0; i < (m_NumDesired+m_NumBurnIn); i++)
   {
      delete [] m_pAccepted[i].x;
   }
   delete [] m_pAccepted;

   for(i = 0; i < m_SamplesPerIter; i++)
   {
      delete [] m_pSamples[i].x;
   }
   delete [] m_pSamples;
   
   m_MaxSamples = 0;
   m_NumDesired = 0;
   m_NumBurnIn = 0;
   m_NumFound = 0;
   m_SamplesPerIter = 0;
   m_MinWSSE = HUGE_VAL;

   delete [] m_pMyBuf;
   delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   delete [] m_pBuf;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
ComputeLikelihood()

Compute the likelihood ratio using the user-selected formualtion.

Stedinger formulation (from Stedinger et al. (2008) “Appraisal of the GLUE Method”): 
   exp[n/2(1-(WSSE/WSSEmin))]
******************************************************************************/
double RejectionSampler::ComputeLikelihoodRatio(double wsse)
{
   double test;
   if(m_bMetropolis == false) // rejection sampling using min. WSSE value
      test = m_MinWSSE;
   else //Metropolis algorithm uses last accepted WSSE value
      test = m_LastWSSE;

   int n = m_pModel->GetObsGroupPtr()->GetNumObs();
   if(m_bStedinger) //formal likelihood ratio
      return exp(((double)n/2.0)*(1.00-(wsse/test)));
   else if (m_bBeven) //pseudo likelihood ratio
      return pow((wsse/test),-1.00*m_ShapeFactor);
   else //default to Stedinger
      return exp(((double)n/2.0)*(1.00-(wsse/test)));
}/* end ComputeLikelihoodRatio() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using a rejection sampler.
******************************************************************************/
void RejectionSampler::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using a rejection sampler. Only useful if
objective function is WSSE.
******************************************************************************/
void RejectionSampler::Optimize(void)
{
   int num,id,i,j,m,g, numSamples, iBest;
   double range, r, rval, upr, lwr;
   int maxGens;
   StatusStruct pStatus;
   ParameterGroup * pGroup;
   
   InitFromFile(GetInFileName());

   maxGens = 1+(int)((m_MaxSamples/(long long)m_SamplesPerIter));

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if(id == 0)
   {
      if(m_bMetropolis == false)
      {
         WriteSetup(m_pModel, "Rejection Sampler");
      }
      else //(bMetropolis == true)
      {
         WriteSetup(m_pModel, "Metropolis MCMC Sampler");
      }

      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Num Found");
   }/* end if() */

   //allocate list of behavioral samples
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();

   NEW_PRINT("SampleStruct", (m_NumDesired+m_NumBurnIn));
   m_pAccepted = new SampleStruct[(m_NumDesired+m_NumBurnIn)];
   MEM_CHECK(m_pAccepted);

   for(i = 0; i < (m_NumDesired+m_NumBurnIn); i++)
   {
      NEW_PRINT("double", num);
      m_pAccepted[i].x = new double[num];
      m_pAccepted[i].fx = HUGE_VAL;
      m_pAccepted[i].n = num;      
   }/* end for() */
   MEM_CHECK(&(m_pAccepted[i-1]));

   //allocate list of random samples
   NEW_PRINT("SampleStruct", m_SamplesPerIter);
   m_pSamples = new SampleStruct[m_SamplesPerIter];
   MEM_CHECK(m_pSamples);

   for(i = 0; i < m_SamplesPerIter; i++)
   {
      NEW_PRINT("double", num);
      m_pSamples[i].x = new double[num];
      m_pSamples[i].fx = HUGE_VAL;
      m_pSamples[i].n = num;      
   }/* end for() */
   MEM_CHECK(&(m_pSamples[i-1]));
 
   //main optimization loop   
   for(g = 0; g < maxGens; g++)
   {
      pStatus.curIter = g+1;
      if(IsQuit() == true){ break;}
      if(m_NumFound >= (m_NumDesired+m_NumBurnIn)){ pStatus.pct = 100.00; break;}
      numSamples = g*m_SamplesPerIter;
      if(numSamples >= m_MaxSamples){ pStatus.pct = 100.00; break;} 

      if(id == 0)
      {
         /*-----------------------------------------------------------
         Generate random samples
         -----------------------------------------------------------*/
         for(i = 0; i < m_SamplesPerIter;  i++)
         {
            for(j = 0; j < num; j++) //for each parameter
            {
               lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
               upr = pGroup->GetParamPtr(j)->GetUprBnd();

               range = upr - lwr;
               r = (double)MyRand() / (double)MY_RAND_MAX;
               rval = (r * range) + lwr;
               m_pSamples[i].x[j] = rval;               
            }/* end for() */
         }/* end for() */
      } /* end if() */

      //evaluate samples, possibly in parallel
      EvaluateSamples();

      //printf("Samples evaluated\n");
   
      if(id == 0)
      {
         /* -----------------------------------------
         apply probabilistic acceptance criteria to
         each sample
         ----------------------------------------- */         
         double r, p;

         for(i = 0; i < m_SamplesPerIter;  i++)
         {
            r = UniformRandom();
            p = ComputeLikelihoodRatio(m_pSamples[i].fx);

            if(r < p)
            {
               //insert new entry
               m_pAccepted[m_NumFound].fx = m_pSamples[i].fx;
               for(m = 0; m < num; m++)
               {
                  m_pAccepted[m_NumFound].x[m] = m_pSamples[i].x[m];
               }
               m_NumFound++;
               m_LastWSSE = m_pSamples[i].fx; //update last accepted WSSE (used by Metropolis MCMC)
               if(m_NumFound >= (m_NumDesired+m_NumBurnIn)) break;
            }/* end if(r < p) */
         }/* end for() */
      }/* end if(id==0) */

      //synchronize num found
      MPI_Bcast(&m_NumFound, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      if(id == 0)
      {
         //determine best entry of the latest sample
         iBest = 0;
         for(j = 1; j < m_SamplesPerIter; j++)
         {
            if(m_pSamples[j].fx < m_pSamples[iBest].fx)
            {                      
               iBest = j;
            }/* end if() */
         }/* end for() */
      
         pGroup->WriteParams(m_pSamples[iBest].x);         

         //constrict parameter bounds
         if((m_TelescopeRate < 1.00) && (m_TelescopeRate > 0.00))
         {
            double pmax, pmin, pval, rate;
            rate = m_TelescopeRate;
            if(m_NumFound > 0) //must have at least one accepted entry
            {
               for(m = 0; m < num; m++) //for each parameter
               {
                  pmax = lwr = pGroup->GetParamPtr(m)->GetLwrBnd();
                  pmin = upr = pGroup->GetParamPtr(m)->GetUprBnd();
                  //determine range of accepted values for the parameter 
                  for(i = 0; i < m_NumFound;  i++) //for each entry
                  {
                     pval = m_pAccepted[i].x[m];
                     if(pval > pmax) pmax = pval;
                     if(pval < pmin) pmin = pval;
                  }/* end for() */
                  lwr += (pmin - lwr)*rate;
                  upr -= (upr - pmax)*rate; 
                  pGroup->GetParamPtr(m)->SetLwrBnd(lwr);
                  pGroup->GetParamPtr(m)->SetUprBnd(upr);               
               }/* end for() */
            }/* end if() */
         }/* end if() */

         pStatus.numRuns = m_pModel->GetCounter();
         pStatus.pct = ((float)100.00*(float)(g+1))/(float)maxGens;
         WriteStatus(&pStatus);
         WriteRecord(m_pModel, (g+1), m_pSamples[iBest].fx, m_NumFound);
      }/* end if(id == 0) */
   }/* end for() */

   if(id == 0)
   { 
      //place model at optimal prameter set
      pGroup->WriteParams(m_pSamples[iBest].x);
      m_pModel->Execute();

      WriteOptimal(m_pModel, m_pSamples[iBest].fx);
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
void RejectionSampler::WriteMetrics(FILE * pFile) 
{
   ParameterGroup * pGroup;
   int nmax, burn, noburn;
   pGroup = m_pModel->GetParamGroupPtr();

   fprintf(pFile, "\nAlgorithm Metrics\n");
   if(m_bMetropolis == false)
   {
      fprintf(pFile, "Algorithm                : Rejection Sampler\n");
   }
   else //if(m_bMetropolis == true)
   {
      fprintf(pFile, "Algorithm                : Metropolis MCMC Sampler\n");
   }

   if(m_bStedinger)
   {
      fprintf(pFile, "Likelihood Function      : Stedinger (formal likelihood)\n");
   }
   else
   {
      fprintf(pFile, "Likelihood Function      : Beven (pseudo likelihood)\n");
      fprintf(pFile, "Shaping Factor           : %0.2lf\n", m_ShapeFactor);
   }
   fprintf(pFile, "Min. WSSR                : %E\n", m_MinWSSE);
   fprintf(pFile, "Max Samples              : %lld\n", m_MaxSamples);
   fprintf(pFile, "Desired Burn In Samples  : %d\n", m_NumBurnIn);
   fprintf(pFile, "Desired Accepted Samples : %d (including burn-in)\n", m_NumDesired + m_NumBurnIn); 
   fprintf(pFile, "Desired Accepted Samples : %d (not including burn-in)\n", m_NumDesired); 
   burn = m_NumFound;
   if(burn > m_NumBurnIn) burn = m_NumBurnIn;
   fprintf(pFile, "Actual Burn In Samples   : %d\n", burn);
   fprintf(pFile, "Actual Accepted Samples  : %d (including burn-in)\n", m_NumFound);
   noburn = m_NumFound - burn;
   if(noburn < 0) noburn = 0;
   fprintf(pFile, "Actual Accepted Samples  : %d (not including burn-in)\n", noburn);
   
   fprintf(pFile,"\nBurn_In_Sample  obj.function  ");
   pGroup->Write(pFile, WRITE_BNR);
   fprintf(pFile,"\n");

   if((m_NumBurnIn <= 0) || (m_NumFound <= 0))
      fprintf(pFile,"no burn in samples were collected\n\n");

   if(m_NumFound > m_NumBurnIn) nmax = m_NumBurnIn;
   else nmax = m_NumFound;
   for(int i = 0; i < nmax; i++)
   {
	   fprintf(pFile, "%-4d            ", i);
      WritePreciseNumber(pFile, m_pAccepted[i].fx);
      fprintf(pFile, "  ");
      pGroup->WriteParams(m_pAccepted[i].x);
      pGroup->Write(pFile, WRITE_SCI);
      fprintf(pFile, "\n");
   }

   fprintf(pFile,"\nAccepted_Sample  obj.function  ");
   pGroup->Write(pFile, WRITE_BNR);
   fprintf(pFile,"\n");

   if(nmax >= m_NumFound)
      fprintf(pFile,"no accepted samples were collected\n");

   for(int i = nmax; i < m_NumFound; i++)
   {
	   fprintf(pFile, "%-4d             ", (i-nmax)+1);
      WritePreciseNumber(pFile, m_pAccepted[i].fx);
      fprintf(pFile, "  ");
      pGroup->WriteParams(m_pAccepted[i].x);
      pGroup->Write(pFile, WRITE_SCI);
      fprintf(pFile, "\n");
   }
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
EvaluateSamples()

Evaluates the objective function of each sample.
******************************************************************************/
void RejectionSampler::EvaluateSamples(void)
{   
   int i, n;   
   ParameterGroup * pGroup;
   double val;

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   
   if(n == 1) //serial execution
   {
      if(IsSuperMUSE() == false)
      {
         WriteInnerEval(WRITE_GLUE, m_SamplesPerIter, '.');
         pGroup = m_pModel->GetParamGroupPtr();
         for(i = 0; i < m_SamplesPerIter; i++) 
         { 
            WriteInnerEval(i+1, m_SamplesPerIter, '.');
            pGroup->WriteParams(m_pSamples[i].x);

            val = m_pModel->Execute();
            m_pSamples[i].fx = val;
         }
         WriteInnerEval(WRITE_ENDED, m_SamplesPerIter, '.');
      }
      else
      {
         EvalSamplesSuperMUSE();
      }
   }/* end if() */
   else /* parallel execution */
   {
      BcastSamples();
      //EvalSamplesParallelAsync();      
      EvalSamplesParallel();
   }/* end else() */
} /* end EvaluateSamples() */

/******************************************************************************
EvalSamplesSuperMUSE()

Compute objective functions of the samples using SuperMUSE. This routine 
interfaces with the RepeatTasker SuperMUSE program, which assigns model 
evaluations to SuperMUSE clients on a first-come-first-served basis.
******************************************************************************/
void RejectionSampler::EvalSamplesSuperMUSE(void)
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
   for(i = 0; i < m_SamplesPerIter; i++)
   {
      //stuff the parameter group with values
      pGroup->WriteParams(m_pSamples[i].x);
         
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
      EvaluateSamples();
   }
   else //SuperMUSE was successful
   {
      for(i = 0; i < m_SamplesPerIter; i++)
      {
         /* -----------------------------------------------
         Stuff the parameter group with ith sample 
         result. This ensures that each objective function 
         gets associated with the correct parameter values.
         ------------------------------------------------ */
         pGroup->WriteParams(m_pSamples[i].x);

         //stuff i-th result into result vector
         val = pSMUSE->GatherResult(i);
         m_pSamples[i].fx = val;
      }/* end for() */
   }/* end else() */
}/* end EvalSamplesSuperMUSE() */

/******************************************************************************
BcastSamples()

When in parallel, only the master computes the samples. All the other 
processors just compute the objeective functions. The BcastSamples() routine is 
called upon to broadcast  the current set of samples from the master processor 
to all of the slave processors.
******************************************************************************/
void RejectionSampler::BcastSamples(void)
{
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix
   pop_size = m_SamplesPerIter;
   num_vars = m_pSamples[0].n;

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
         m_pBuf[idx] = m_pSamples[j].x[i];
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
         m_pSamples[j].x[i] = m_pBuf[idx];
      }/* end for() */
   }/* end for() */
}/* end BcastSamples() */

/******************************************************************************
EvalSamplesParallel()

Compute objective function of entire set of samples in parallel. Each processor 
evaluates a predetermined number of samples, based on their processor id.
******************************************************************************/
void RejectionSampler::EvalSamplesParallel(void)
{    
   int i ,j, num_procs, id, bufsize, idx, proc;
   ParameterGroup * pGroup;
   MPI_Status status;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
   bufsize = (m_SamplesPerIter/num_procs) + 1;

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];

      NEW_PRINT("double", bufsize);
      m_pTmpBuf = new double[bufsize];

      NEW_PRINT("double", m_SamplesPerIter);
      m_pBigBuf = new double[m_SamplesPerIter];
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_SamplesPerIter; i++) 
   { 
      if((i % num_procs) == id)
      {          
         pGroup->WriteParams(m_pSamples[i].x);

         m_pMyBuf[j] = m_pModel->Execute();
         m_pTmpBuf[j] = m_pMyBuf[j];
         j++;
      }/* end if() */
   }/* end for() */

   //gather results
#if 0
   // old and inefficient way used MPI_Bcast to effect the gather operation!
   for(i = 0; i < num_procs; i++)
   {
	   //receive someones buf, this will clobber myBuf
      MPI_Bcast(m_pMyBuf, bufsize, MPI_DOUBLE, i, MPI_COMM_WORLD);

      for(j = 0; j < bufsize; j++)
      {
         idx = (num_procs * j) + i; //idx maps myBuf into bigBuf

         if(idx < m_SamplesPerIter)
         {
            m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            m_pMyBuf[j] = m_pTmpBuf[j]; //restore myBuf...clobbered by bcast
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //stuff results into swarm
   for(i = 0; i < m_SamplesPerIter; i++)
   {
      m_pSamples[i].fx = m_pBigBuf[i];
   }/* end for() */
#endif
   //subordinates - send results to supervisor
   if(id != 0)
   {
      MPI_Send(m_pMyBuf, bufsize, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);
   }/* end if(subordinate) */
   //supervisor - receive subordinate results
   else
   {
      for (proc = 1; proc < num_procs; proc++) 
      {
         MPI_Recv(m_pMyBuf, bufsize, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_RESULTS_TAG, MPI_COMM_WORLD, &status);
         i = status.MPI_SOURCE;

         for(j = 0; j < bufsize; j++)
         {
            idx = (num_procs * j) + i; //idx maps myBuf into bigBuf

            if(idx < m_SamplesPerIter)
            {
               m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            }/* end if() */
         }/* end for() */
      } /* end for(each subordinate) */

      //interleave supervisor results
      for(j = 0; j < bufsize; j++)
      {
         idx = (num_procs * j); //idx maps myBuf into bigBuf

         if(idx < m_SamplesPerIter)
         {
            m_pBigBuf[idx] = m_pTmpBuf[j]; //gather into bigbuf
         }/* end if() */
      }/* end for() */

      //stuff results into swarm
      for(i = 0; i < m_SamplesPerIter; i++)
      {
         m_pSamples[i].fx = m_pBigBuf[i];
      }/* end for() */
   }/* end else(supervisor) */

   MPI_Barrier(MPI_COMM_WORLD);

}/* end EvalSamplesParallel() */

/******************************************************************************
EvalSamplesParallelAsync()

Compute objective function of entire set of samples in parallel. Compute tasks 
are distributed asynchronously by the supervisor processor.
******************************************************************************/
void RejectionSampler::EvalSamplesParallelAsync(void)
{    
   int i ,j, num_procs, num_subords, id, bufsize, proc;
   ParameterGroup * pGroup;
   MPI_Status status;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   num_subords = num_procs - 1;
   
   bufsize = m_SamplesPerIter;

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];
   }

   //perform parallel evaluations
   pGroup = m_pModel->GetParamGroupPtr();
   int nleft = m_SamplesPerIter;
   int ntasks = 0;
   double fx;
   while(1)
   {
      //supervisor waits for task requests from subordinates
      if(id == 0)
      {
         MPI_Recv(&proc, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_REQUEST_TAG, MPI_COMM_WORLD, &status);
         i = m_SamplesPerIter - nleft;
         MPI_Send(&i, 1, MPI_INTEGER, proc, MPI_INDEX_TAG, MPI_COMM_WORLD);
         nleft--;
         /* --------------------------------------------------------
         when this if statement is satisfied, all subordinates will 
         have been indirectly informed that there are no more tasks
         ----------------------------------------------------------- */
         if(nleft <= -1*num_subords)
         {
            break;
         }/* end if() */
      }/* end if( supervisor ) */
      //subordinates request work from supervisor
      else
      {
         MPI_Send(&id, 1, MPI_INTEGER, 0, MPI_REQUEST_TAG, MPI_COMM_WORLD);
         MPI_Recv(&i, 1, MPI_INTEGER, 0, MPI_INDEX_TAG, MPI_COMM_WORLD, &status);
         //if not indexable, then there are no more tasks
         if(i >= m_SamplesPerIter)
         {
            break;
         }
         //process the task
         else
         {
            m_pMyBuf[ntasks] = (double)i;
            pGroup->WriteParams(m_pSamples[i].x);
            m_pSamples[i].fx = m_pModel->Execute();
            ntasks++;
         }/* end else() */
      }/* end else ( subordinates ) */
   }/* end while() */

   //printf("Processor %d is ready to gather results\n", id);

   //gather results
   if(id == 0)
   {
      //collect results from suboordinates
      for(j = 0; j < m_SamplesPerIter; j++)
      {
         MPI_Recv(&i, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_INDEX_TAG, MPI_COMM_WORLD, &status);
         proc = status.MPI_SOURCE;
         MPI_Recv(&fx, 1, MPI_DOUBLE, proc, MPI_RESULTS_TAG, MPI_COMM_WORLD, &status);
         m_pSamples[i].fx = fx;
      }/* end for() */
   }/* end if() */
   else
   {
      for(j = 0; j < ntasks; j++)
      {
         i = (int)m_pMyBuf[j];
         fx = m_pSamples[i].fx;
         MPI_Send(&i, 1, MPI_INTEGER, 0, MPI_INDEX_TAG, MPI_COMM_WORLD);
         MPI_Send(&fx, 1, MPI_DOUBLE, 0, MPI_RESULTS_TAG, MPI_COMM_WORLD);
      }/* end for() */
   }/* end else() */

   //printf("\nProcessor %d is waiting at the barrier\n", id);
   MPI_Barrier(MPI_COMM_WORLD);
}/* end EvalSamplesParallelAsync() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void RejectionSampler::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ], type[DEF_STR_SZ];

   char startToken[DEF_STR_SZ];
   char endToken[DEF_STR_SZ];

   m_MaxSamples = 100;
   m_NumDesired = 10;
   m_NumFound = 0;
   m_SamplesPerIter = 10;
   m_MinWSSE = HUGE_VAL;
   m_bStedinger = true;
   m_bBeven = false;
   m_ShapeFactor = 0.5; //RMSE

   if(m_bMetropolis == false)
   {
      strcpy(startToken, "BeginRejectionSampler");
      strcpy(endToken, "EndRejectionSampler");
   }
   else //if(m_bMetropolis == true)
   {
      strcpy(startToken, "BeginMetropolisSampler");
      strcpy(endToken, "EndMetropolisSampler");
   }

   //read in rejection sampling configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, startToken, pFileName) == true)
   {
      FindToken(pFile, endToken, pFileName);
      rewind(pFile);

      FindToken(pFile, startToken, pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, endToken) == NULL)
      {         
         if(strstr(line, "SamplesPerIter") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_SamplesPerIter); 
            if(m_SamplesPerIter < 1)
            {
               LogError(ERR_FILE_IO, "Invalid setting. Defaulting to 10.");
               m_SamplesPerIter = 10;
            }
         }/*end if() */         
         else if(strstr(line, "NumDesired") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumDesired); 
            if(m_NumDesired < 1)
            {
               LogError(ERR_FILE_IO, "Invalid setting. Defaulting to 10.");
               m_NumDesired = 10;
            }
         }/*end else if() */         
         else if(strstr(line, "BurnInSamples") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumBurnIn); 
            if(m_NumBurnIn < 0)
            {
               LogError(ERR_FILE_IO, "Invalid setting. Defaulting to 0 (no burn in).");
               m_NumBurnIn = 0;
            }
         }/*end else if() */
         else if(strstr(line, "MaxSamples") != NULL)
         {
            sscanf(line, "%s %lld", tmp, &m_MaxSamples); 
            if(m_MaxSamples < 1)
            {
               LogError(ERR_FILE_IO, "Invalid setting. Defaulting to 100.");
               m_MaxSamples = 100;
            }
         }/*end else if() */         
         else if(strstr(line, "MinWSSE") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_MinWSSE); 
         }/*end else if() */
         else if(strstr(line, "LikelihoodType") != NULL)
         {
            sscanf(line, "%s %s", tmp, type);
			   MyStrLwr(type);
			   if(strcmp(type, "stedinger") == 0)
            {
               m_bStedinger = true;
               m_bBeven = false;
            }
			   else if(strcmp(type, "beven") == 0)
            {
               m_bStedinger = false;
               m_bBeven = true;
            }
			   else 
			   {
			      sprintf(tmp, "Unknown Likelihood Type |%s|, defaulting to Stedinger", type);
               LogError(ERR_FILE_IO, tmp);
			   }
         }/*end else if() */
         else if(strstr(line, "ShapingFactor") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_ShapeFactor);
         }/*end else if() */
         else if(strstr(line, "TelescopeRate") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_TelescopeRate);
         }/*end else if() */
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   
   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
RJSMP_Program()

Calibrate or optimize the model using Rejection Sampling.
******************************************************************************/
void RJSMP_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("RejectionSampler", 1);
   RejectionSampler * RJSMP = new RejectionSampler(model, false);
   MEM_CHECK(RJSMP);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { RJSMP->Calibrate(); }
   else { printf("Rejection Sampling algorithm can only be used with WSSE objective function.\n"); }

   delete RJSMP;
   delete model;
} /* end RJSMP_Program() */

/******************************************************************************
METRO_Program()

Calibrate or optimize the model using Metroplois MCMS Sampling.
******************************************************************************/
void METRO_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("RejectionSampler", 1);
   RejectionSampler * METRO = new RejectionSampler(model, true);
   MEM_CHECK(METRO);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { METRO->Calibrate(); }
   else { printf("Metropolis MCMC Sampling algorithm can only be used with WSSE objective function.\n"); }

   delete METRO;
   delete model;
} /* end METRO_Program() */

