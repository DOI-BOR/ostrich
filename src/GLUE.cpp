/******************************************************************************
File     : GLUE.cpp
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Generalized Likelihood Uncertainty Engine - GLUE

Version History
06-23-10    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "GLUE.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "Utility.h"

double mpi_timer_count = 0;
double mpi_timer_start = 0;
double mpi_timer_end = 0;
double gTotal = 0;
double eTotal = 0;
double uTotal = 0;
double oTime = 0;
double iTime = 0;

int gSerialCount = 0;

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
GLUE::GLUE(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pBehavioral = NULL;
   m_pSamples = NULL;
   m_MaxSamples = 0;
   m_NumDesired = 0;
   m_NumFound = 0;
   m_SamplesPerIter = 0;
   m_Threshold = -1.00;

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   //m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   //m_pBuf    = NULL;
   m_iCounts = NULL;
   m_iDispls = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the GLUE and it's member variables.
******************************************************************************/
void GLUE::Destroy(void)
{
   int i;
   for(i = 0; i < m_NumDesired; i++)
   {
      delete [] m_pBehavioral[i].x;
   }
   delete [] m_pBehavioral;

   for(i = 0; i < m_SamplesPerIter; i++)
   {
      delete [] m_pSamples[i].x;
   }
   delete [] m_pSamples;
   
   m_MaxSamples = 0;
   m_NumDesired = 0;
   m_NumFound = 0;
   m_SamplesPerIter = 0;
   m_Threshold = -1.00;

   delete [] m_pMyBuf;
   //delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   //delete [] m_pBuf;
   delete [] m_iCounts;
   delete [] m_iDispls;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using GLUE.
******************************************************************************/
void GLUE::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using GLUE.
******************************************************************************/
void GLUE::Optimize(void)
{
   double oStart, oEnd;
   double tStart, tEnd;
   int num,id,i,j,k,m,g,nprocs;
   long long numSamples;
   double range, r, rval, upr, lwr;
   int maxGens;
   StatusStruct pStatus;
   ParameterGroup * pGroup;
   
   oStart = GetElapsedTics();

   InitFromFile(GetInFileName());

   maxGens = 1+(int)(m_MaxSamples/(long long)m_SamplesPerIter);

   mpi_timer_start = GetElapsedTics();
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   double divisor = (double) m_SamplesPerIter / (double)nprocs;
   m_iStart=(int)(ceil(divisor*(double)id));
   m_iEnd=(int)(ceil(divisor*(double)(id+1)));
   if(m_iEnd > m_SamplesPerIter) m_iEnd = m_SamplesPerIter;

   //collect size and offset of each processors evaluation buffer, for synhronizing subsequent gather operations
   int iSize = (m_iEnd - m_iStart);
   m_iCounts = new int[nprocs];
   m_iDispls = new int[nprocs];
   MPI_Gather(&iSize, 1, MPI_INTEGER, m_iCounts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
   MPI_Gather(&m_iStart, 1, MPI_INTEGER, m_iDispls, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

   mpi_timer_end = GetElapsedTics();
   mpi_timer_count += (mpi_timer_end - mpi_timer_start);
    
   if(id == 0)
   {
      WriteSetup(m_pModel, "Generalized Likelihood Uncertainty Engine");
      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Num Found");
   }/* end if() */

   //allocate list of behavioral samples
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();

   NEW_PRINT("SampleStruct", m_NumDesired);
   m_pBehavioral = new SampleStruct[m_NumDesired];
   MEM_CHECK(m_pBehavioral);

   for(i = 0; i < m_NumDesired; i++)
   {
      NEW_PRINT("double", num);
      m_pBehavioral[i].x = new double[num];
      m_pBehavioral[i].fx = HUGE_VAL;
      m_pBehavioral[i].n = num;      
   }/* end for() */
   MEM_CHECK(&(m_pBehavioral[i-1]));

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
      for(j = 0; j < num; j++) m_pSamples[i].x[j] = -999.999;
   }/* end for() */
   MEM_CHECK(&(m_pSamples[i-1]));
 
   iTime = GetElapsedTics() - oStart;

   //main optimization loop   
   for(g = 0; g < maxGens; g++)
   {
      pStatus.curIter = g+1;
      if(IsQuit() == true){ break;}
      if(m_NumFound == m_NumDesired){ pStatus.pct = 100.00; break;}
      numSamples = (long long)g*(long long)m_SamplesPerIter;
      if(numSamples >= m_MaxSamples){ pStatus.pct = 100.00; break;} 

      tStart = GetElapsedTics();
      /*-----------------------------------------------------------
      Generate random samples
      -----------------------------------------------------------*/
      for(i = m_iStart; i < m_iEnd;  i++)
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
      tEnd = GetElapsedTics();
      gTotal += (tEnd - tStart);

      //evaluate samples, possibly in parallel
      tStart = GetElapsedTics();
      EvaluateSamples();
      tEnd = GetElapsedTics();
      eTotal += tEnd - tStart;
  
      //revise list of best samples (master only)
      tStart = GetElapsedTics();
      if(id == 0)
      {
         for(i = 0; i < m_SamplesPerIter; i++)
         {
            for(j = 0; j < m_NumDesired; j++)
            {
               if(m_pSamples[i].fx < m_pBehavioral[j].fx)
               {
                  //shift previous list to make room
                  for(k = m_NumDesired-1; k > j; k--)
                  {
                     for(m = 0; m < num; m++)
                     {
                        m_pBehavioral[k].x[m] = m_pBehavioral[k-1].x[m];
                        m_pBehavioral[k].fx = m_pBehavioral[k-1].fx;
                     }
                  }
                  //insert new entry
                  for(m = 0; m < num; m++)
                  {
                     m_pBehavioral[j].x[m] = m_pSamples[i].x[m];
                     m_pBehavioral[j].fx = m_pSamples[i].fx;
                  }
                  break;
               }/* end if() */
            }/* end for() */
         }/* end for() */      

         //count number of behavioral
         m_NumFound = 0;
         for(j = 0; j < m_NumDesired; j++)
         {
            if(m_pBehavioral[j].fx < m_Threshold)
            {                      
               m_NumFound++;
            }/* end if() */
         }/* end for() */
      }/* end if() */
      mpi_timer_start = GetElapsedTics();
      MPI_Bcast(&m_NumFound, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      tEnd = mpi_timer_end = GetElapsedTics();
      mpi_timer_count += (tEnd - mpi_timer_start);
      uTotal += (tEnd - tStart);

      //first entry is always the best
      if(id == 0) { pGroup->WriteParams(m_pBehavioral[0].x); }

      pStatus.pct = ((float)100.00*(float)(g+1))/(float)maxGens;
      pStatus.numRuns = m_pModel->GetCounter();
      if(id == 0){ WriteStatus(&pStatus); }
      if(id == 0){ WriteRecord(m_pModel, (g+1), m_pBehavioral[0].fx, m_NumFound);}
   }/* end for() */

   if(id == 0)
   { 
      //place model at optimal prameter set
      pGroup->WriteParams(m_pBehavioral[0].x);
      m_pModel->Execute();

      WriteOptimal(m_pModel, m_pBehavioral[0].fx);
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      //write algorithm metrics
      //WriteAlgMetrics(this);
   }

   oEnd = GetElapsedTics();
   oTime = (oEnd - oStart);

   if(id == 0)
   { 
      WriteAlgMetrics(this);
   }
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void GLUE::WriteMetrics(FILE * pFile) 
{
   ParameterGroup * pGroup;
   pGroup = m_pModel->GetParamGroupPtr();

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Generalized Likelihood Uncertainty Estimation\n");
   fprintf(pFile, "Behavioral Threshold    : %E\n", m_Threshold);
   fprintf(pFile, "Max Samples             : %lld\n", m_MaxSamples);
   fprintf(pFile, "Actual Num. Behavorial  : %d\n", m_NumFound);
   fprintf(pFile, "Desired Num. Behavorial : %d\n", m_NumDesired);
   fprintf(pFile, "Seconds Spent on MPI    : %lf\n", mpi_timer_count);
   fprintf(pFile, "Secs Spent in Setup     : %lf\n", iTime);
   fprintf(pFile, "Secs Spent in Evaluate  : %lf\n", eTotal);
   fprintf(pFile, "Secs Spent in Generate  : %lf\n", gTotal);
   fprintf(pFile, "Secs Spent in Update    : %lf\n", uTotal);
   fprintf(pFile, "Secs required overall   : %lf\n", oTime);
   fprintf(pFile, "Secs not accounted for  : %lf\n\n", oTime - eTotal - gTotal - uTotal - iTime);
   fprintf(pFile, "Serial Count            : %d\n\n", gSerialCount);

   fprintf(pFile,"Sample  obj.function  ");
   pGroup->Write(pFile, WRITE_BNR);
   fprintf(pFile,"\n");

   for(int i = 0; i < m_NumDesired; i++)
   {
	   fprintf(pFile, "%-4d  ", i);
      WritePreciseNumber(pFile, m_pBehavioral[i].fx);
      fprintf(pFile, "  ");
      pGroup->WriteParams(m_pBehavioral[i].x);
      pGroup->Write(pFile, WRITE_SCI);
      fprintf(pFile, "\n");
   }
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
EvaluateSamples()

Evaluates the objective function of each sample.
******************************************************************************/
void GLUE::EvaluateSamples(void)
{   
   int i, n;   
   ParameterGroup * pGroup;
   double val;

   mpi_timer_start = GetElapsedTics();
   MPI_Comm_size(MPI_COMM_WORLD, &n);
   mpi_timer_end = GetElapsedTics();
   mpi_timer_count += (mpi_timer_end - mpi_timer_start);
   
   if(n == 1) //serial execution
   {
      gSerialCount++;

      //WriteInnerEval(WRITE_GLUE, m_SamplesPerIter, '.');
      pGroup = m_pModel->GetParamGroupPtr();
      for(i = 0; i < m_SamplesPerIter; i++) 
      { 
         //WriteInnerEval(i+1, m_SamplesPerIter, '.');
         pGroup->WriteParams(m_pSamples[i].x);

         val = m_pModel->Execute();
         m_pSamples[i].fx = val;
      }
      //WriteInnerEval(WRITE_ENDED, m_SamplesPerIter, '.');
   }/* end if() */
   else /* parallel execution */
   {
      //BcastSamples();
      EvalSamplesParallel();      
   }/* end else() */
} /* end EvaluateSamples() */

/******************************************************************************
BcastSamples()

When in parallel, only the master computes the samples. All the other 
processors just compute the objeective functions. The BcastSamples() routine is 
called upon to broadcast  the current set of samples from the master processor 
to all of the slave processors.
******************************************************************************/
#if 0
void GLUE::BcastSamples(void)
{
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   mpi_timer_start = GetElapsedTics();
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   mpi_timer_end = GetElapsedTics();
   mpi_timer_count += (mpi_timer_end - mpi_timer_start);

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
   mpi_timer_start = GetElapsedTics();
   MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   mpi_timer_end = GetElapsedTics();
   mpi_timer_count += (mpi_timer_end - mpi_timer_start);

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
#endif

/******************************************************************************
EvalSamplesParallel()

Compute objective function of entire set of samples in parallel. Each processor 
evaluates a predetermined number of samples, based on their processor id.
******************************************************************************/
void GLUE::EvalSamplesParallel(void)
{    
   int i ,j, k, num_procs, id, bufsize;
   ParameterGroup * pGroup;

   //setup processor id and number of processors
   mpi_timer_start = GetElapsedTics();
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   mpi_timer_end = GetElapsedTics();
   mpi_timer_count += (mpi_timer_end - mpi_timer_start);
   
   bufsize = (m_iEnd - m_iStart);

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];

      //NEW_PRINT("double", bufsize);
      //m_pTmpBuf = new double[bufsize];

      NEW_PRINT("double", m_SamplesPerIter);
      m_pBigBuf = new double[m_SamplesPerIter];
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = m_iStart; i < m_iEnd; i++) 
   { 
      pGroup->WriteParams(m_pSamples[i].x);
      m_pMyBuf[j] = m_pModel->Execute();
      j++;
   }/* end for() */

   //gather F(x) results at master
   mpi_timer_start = GetElapsedTics();
   MPI_Gatherv(m_pMyBuf, bufsize, MPI_DOUBLE, m_pBigBuf, m_iCounts, m_iDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   mpi_timer_end = GetElapsedTics();
   mpi_timer_count += (mpi_timer_end - mpi_timer_start);

   //stuff results into sample list
   for(i = 0; i < m_SamplesPerIter; i++)
   {
      m_pSamples[i].fx = m_pBigBuf[i];
   }/* end for() */

   for(j = 0; j < m_pSamples[0].n; j++)
   {
      k = 0;
      for(i = m_iStart; i < m_iEnd; i++)
      {
         m_pMyBuf[k] = m_pSamples[i].x[j];
         k++;
      }/* end for() */

      //gather parameter values at master
      mpi_timer_start = GetElapsedTics();
      MPI_Gatherv(m_pMyBuf, bufsize, MPI_DOUBLE, m_pBigBuf, m_iCounts, m_iDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      mpi_timer_end = GetElapsedTics();
      mpi_timer_count += (mpi_timer_end - mpi_timer_start);

      //stuff results into sample list
      for(i = 0; i < m_SamplesPerIter; i++)
      {
         m_pSamples[i].x[j] = m_pBigBuf[i];
      }/* end for() */
   }/* end for() */
}/* end EvalSamplesParallel() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void GLUE::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];

   m_MaxSamples = 100;
   m_NumDesired = 10;
   m_NumFound = 0;
   m_SamplesPerIter = 10;
   m_Threshold = 1000.00;

   //read in GLUE configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open GLUE config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginGLUE", pFileName) == true)
   {
      FindToken(pFile, "EndGLUE", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginGLUE", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndGLUE") == NULL)
      {         
         if(strstr(line, "SamplesPerIter") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_SamplesPerIter); 
            if(m_SamplesPerIter < 1)
            {
               LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 10.");
               m_SamplesPerIter = 10;
            }
         }/*end if() */         
         else if(strstr(line, "NumBehavioral") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumDesired); 
            if(m_NumDesired < 1)
            {
               LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 10.");
               m_NumDesired = 10;
            }
         }/*end else if() */         
         else if(strstr(line, "MaxSamples") != NULL)
         {            
            sscanf(line, "%s %lld", tmp, &m_MaxSamples); 
            if(m_NumDesired < 1)
            {
               LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 100.");
               m_MaxSamples = 100;
            }
         }/*end else if() */         
         else if(strstr(line, "Threshold") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Threshold); 
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

   SetObjFuncThreshold(m_Threshold);
} /* end InitFromFile() */

/******************************************************************************
GLUE_Program()

Calibrate or optimize the model using GLUE.
******************************************************************************/
void GLUE_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("GLUE", 1);
   GLUE * pGLUE = new GLUE(model);
   MEM_CHECK(pGLUE);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { pGLUE->Calibrate(); }
   else { pGLUE->Optimize(); }

   delete pGLUE;
   delete model;
} /* end GLUE_Program() */

