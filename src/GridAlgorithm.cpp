/******************************************************************************
File     : GridAlgorithm.cpp
Author   : L. Shawn Matott
Copyright: 2005, L. Shawn Matott

This algorithm will "grid" the objective function surface by evaluating an
exhaustive (or nearly exhaustive) set of parameter configurations.

Version History
01-25-05    lsm   created
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel.
07-18-07    lsm   Added support for SuperMUSE
******************************************************************************/
#include "mpi_stub.h"
#include <stdlib.h>
#include <string.h>

#include "GridAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "SuperMUSE.h"
#include "StatsClass.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"
#include "SuperMuseUtility.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void GridAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
GridAlgorithm::GridAlgorithm(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;   
   m_Lwr = NULL;
   m_Rval = NULL;
   m_pMini = NULL;
   m_pStats = NULL;
   m_pDims = NULL;

   m_GridSize = 0;
   m_NumIters = 0;
   m_CurIter = 0;   
   m_MiniSize = 0;
   m_pBest = NULL;

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   m_pBuf    = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the grid algorithm and it's member variables.
******************************************************************************/
void GridAlgorithm::Destroy(void)
{
   int i;

   delete [] m_pDims;
   delete [] m_Lwr;
   delete [] m_Rval;

   //free up mini-grid
   delete [] m_pMini->dp;
   delete [] m_pMini->f;
   for(i = 0; i < m_MiniSize; i++){delete [] (m_pMini->p[i]);}
   delete [] m_pMini->p;
   delete m_pMini;
   
   delete m_pStats;

   m_GridSize = 0;
   m_NumIters = 0;
   m_CurIter = 0;
   m_MiniSize = 0;   
   delete [] m_pBest;

   delete [] m_pMyBuf;
   delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   delete [] m_pBuf;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using Grid algorithm (exhaustive 
search).
******************************************************************************/
void GridAlgorithm::Calibrate(void)
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

Minimize the objective function using Grid algorithm (exhaustive search).
******************************************************************************/
void GridAlgorithm::Optimize(void)
{
   StatusStruct pStatus;
   int num, g, offset, offsize;
   int id, idx;
   int i, j;   
   double median;
   
   InitFromFile(GetInFileName());

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //initialize the mini-grid
   num = m_pMini->nprm;
   m_NumLeft = m_GridSize;
   if(m_NumLeft < m_MiniSize){ m_MiniSize = m_NumLeft; m_NumIters = 1;}
   for(i = 0; i < m_MiniSize; i++) //for each member of the mini-grid
   {
      m_pMini->f[i] = NEARLY_HUGE;
      //for each parameter
      for(j = 0; j < num; j++) {m_pMini->p[i][j] = GetGridVal(i,j); }
   }/* end for() */

   //handle warm start
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
      for(j = 0; j < num; j++) {m_pMini->p[0][j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal(); }      
   }
   //handle parameter extraction
   if(m_pModel->GetParamGroupPtr()->CheckExtraction() == true)
   {
      for(j = 0; j < num; j++) {m_pMini->p[0][j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal(); }
   }

   //write out setup
   if(id == 0)
   {
      WriteSetup(m_pModel, "Grid Algorithm (Exhaustive Search)");
      //write banner
      WriteBanner(m_pModel, "iter   best value     ", "Median Value");
   }/* end if() */

   //main optimization loop   
   pStatus.maxIter = m_NumIters;
   offsize = m_MiniSize;
   for(g = 0; g < m_NumIters; g++)
   {      
      if(IsQuit() == true){ break;}
      pStatus.curIter = m_CurIter = g+1;

      //evaluate (mini) grid, possibly in parallel
      EvaluateGrid();
      m_NumLeft -= m_MiniSize;

      //update grid output with results from mini-grid
      if(id == 0) WriteGrid(m_pMini, m_MiniSize);

      //revise global best
      for(i = 0; i < m_MiniSize; i++)
      {
         if(m_pBest[num] > m_pMini->f[i])
         {
            for(j = 0; j < num; j++) m_pBest[j] = m_pMini->p[i][j];
            m_pBest[j] = m_pMini->f[i];
         }/* end if() */
      }/* end if() */

      median = CalcMedian(m_pMini->f, m_MiniSize);
      m_pModel->GetParamGroupPtr()->WriteParams(m_pBest);

      //update mini-grid with next set of grid values
      offset = (g+1)*offsize;
      if(m_NumLeft < m_MiniSize){ m_MiniSize = m_NumLeft;}
      for(i = 0; i < m_MiniSize; i++) //for each member of the mini-grid
      {
         idx = i+offset;
         m_pMini->f[i] = NEARLY_HUGE;

         //for each parameter
         for(j = 0; j < num; j++) {m_pMini->p[i][j] = GetGridVal(idx,j); }
      }/* end for() */

      if(id == 0)
      { 
         WriteRecord(m_pModel, (g+1), m_pBest[num], median);
         pStatus.pct  = ((float)100.00*(float)m_CurIter)/(float)m_NumIters;
         pStatus.numRuns = m_pModel->GetCounter();
         WriteStatus(&pStatus);
      }

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */

   //place model at optimal prameter set
   m_pModel->GetParamGroupPtr()->WriteParams(m_pBest);
   m_pModel->Execute();

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   { 
      m_MiniSize = offsize;
      WriteOptimal(m_pModel, m_pBest[num]);
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
void GridAlgorithm::WriteMetrics(FILE * pFile)
{
   int i, num_procs;

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Grid Algorithm (Exhaustive Search)\n");
   fprintf(pFile, "Max Iterations          : %d\n", m_NumIters);
   fprintf(pFile, "Actual Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Grid Size               : %d\n", m_GridSize);

   fprintf(pFile, "Grid Dimensions         : ");
   for(i = 0; i < m_pMini->nprm; i++)
   {
      fprintf(pFile, "%d", m_pDims[i]);
      if(i < (m_pMini->nprm - 1)) fprintf(pFile, " by ");
   }
   fprintf(pFile, "\n");

   fprintf(pFile, "Mini Grid Size          : %d\n", m_MiniSize);
   
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   //fprintf(pFile, "Total Evals             : ");
   
   if(num_procs > 1)
   { 
      fprintf(pFile, "Total Evals             : %d\n", m_GridSize - m_NumLeft);
   }
   else
   {  
      m_pModel->WriteMetrics(pFile);
      //fprintf(pFile, "%d\n", m_pModel->GetCounter());
   }
}/* end WriteMetrics() */

/******************************************************************************
EvaluateGrid()

Evaluates the objective function of each set in the grid.
******************************************************************************/
void GridAlgorithm::EvaluateGrid(void)
{   
   int i, n;   
   ParameterGroup * pGroup;
   double val;

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   
   if(n == 1) //serial execution
   {
      if(IsSuperMUSE() == false)
      {
         WriteInnerEval(WRITE_GRID, m_MiniSize, '.');
         pGroup = m_pModel->GetParamGroupPtr();
         for(i = 0; i < m_MiniSize; i++) 
         { 
            WriteInnerEval(i+1, m_MiniSize, '.');
            pGroup->WriteParams(m_pMini->p[i]);
            val = m_pModel->Execute();
            m_pMini->f[i] = val;
         }
         WriteInnerEval(WRITE_ENDED, m_MiniSize, '.');
      }
      else //SuperMUSE
      {
         EvalGridSuperMUSE();
      }
   }/* end if() */
   else /* parallel execution */
   {
      BcastGrid();
      EvalGridParallel();      
   }/* end else() */
} /* end EvaluateGrid() */

/******************************************************************************
BcastGrid()

When in parallel, all processors compute the objeective functions. The 
BcastGrid() routine is called upon to broadcast the current mini grid from the 
master processor to all of the slave processors.
******************************************************************************/
void GridAlgorithm::BcastGrid(void)
{
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix
   pop_size = m_MiniSize;
   num_vars = m_pMini->nprm;

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
         m_pBuf[idx] = m_pMini->p[j][i];
      }/* end for() */
   }/* end for() */

   //broadcast the flattened matrix
   MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //use the flattened matrix to fill the mini grid
   for(i = 0; i < num_vars; i++)
   {
      for(j = 0; j < pop_size; j++)
      {         
         idx = (num_vars)*j + i;
         m_pMini->p[j][i] = m_pBuf[idx];
      }/* end for() */
   }/* end for() */
}/* end BcastGrid() */

/******************************************************************************
EvalGridParallel()

Compute objective function of entire mini-grid in parallel. Each processor 
evaluates a predetermined number of grid members, based on their processor id.
******************************************************************************/
void GridAlgorithm::EvalGridParallel(void)
{    
   int i ,j, num_procs, id, bufsize, idx;
   ParameterGroup * pGroup;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
   bufsize = (m_MiniSize/num_procs) + 1;

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];

      NEW_PRINT("double", bufsize);
      m_pTmpBuf = new double[bufsize];

      NEW_PRINT("double", m_MiniSize);
      m_pBigBuf = new double[m_MiniSize];
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_MiniSize; i++) 
   { 
      if((i % num_procs) == id)
      {          
         pGroup->WriteParams(m_pMini->p[i]);
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

         if(idx < m_MiniSize)
         {
            m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            m_pMyBuf[j] = m_pTmpBuf[j]; //restore myBuf...clobbered by bcast
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //stuff results into mini grid
   for(i = 0; i < m_MiniSize; i++)
   {
      m_pMini->f[i] = m_pBigBuf[i];
   }/* end for() */
}/* end EvalGridParallel() */

/******************************************************************************
EvalGridSuperMUSE()

Compute objective functions of a mini-grid using SuperMUSE. This routine 
interfaces with the RepeatTasker SuperMUSE program, which assigns model 
evaluations to SuperMUSE clients on a first-come-first-served basis.
******************************************************************************/
void GridAlgorithm::EvalGridSuperMUSE(void)
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
   mini-grid member as entries in the task file.

   Entries are first accumlated into a temp file to prevent the 
   SuperMUSE RepeatTasker program from prematurely processing the task
   file.
   ---------------------------------------------------------------- */
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_MiniSize; i++)
   {
      //stuff the parameter group with values
      pGroup->WriteParams(m_pMini->p[i]);
         
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
      EvaluateGrid();
   }
   else //SuperMUSE was successful
   {
      for(i = 0; i < m_MiniSize; i++)
      {
         /* -----------------------------------------------
         Stuff the parameter group with ith grid
         member. This ensures that each objective function 
         gets associated with the correct parameter values.
         ------------------------------------------------ */
         pGroup->WriteParams(m_pMini->p[i]);

         //stuff i-th result into chromosome pool
         val = pSMUSE->GatherResult(i);
         m_pMini->f[i] = val;
      }/* end for() */
   }/* end else() */
}/* end EvalGridSuperMUSE() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void GridAlgorithm::InitFromFile(IroncladString pFileName)
{
   const char * start_tag = "BeginGridAlg";
   const char * end_tag   = "EndGridAlg";
   int    start_len = (int)strlen(start_tag);
   int    end_len   = (int)strlen(end_tag);

   FILE * pFile;
   int i, j, k, num, num_procs;
   double dim, upr, lwr;
   char * pTok;
   char * line;   
   char tmp[DEF_STR_SZ];
   
   num = m_pModel->GetParamGroupPtr()->GetNumParams();

   //set default mini-grid size
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   m_MiniSize = 10*num*num_procs;   

   //allocate space for grid dimensions
   NEW_PRINT("int", num);
   m_pDims = new int[num];
   MEM_CHECK(m_pDims);

   //allocate space for lower bounds
   NEW_PRINT("double", num);
   m_Lwr = new double[num];
   MEM_CHECK(m_Lwr);

   //allocate space for best
   NEW_PRINT("double", num+1);
   m_pBest = new double[num+1];
   MEM_CHECK(m_pBest);
   m_pBest[num] = NEARLY_HUGE; //last index stores f(p)

   //assign default grid dimensions
   for(k = 0; k < num; k++) m_pDims[k] = 2;

   //read in Grid configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open Grid config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, start_tag, pFileName) == true)
   {
      FindToken(pFile, end_tag, pFileName);
      rewind(pFile);

      FindToken(pFile, start_tag, pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strncmp(line, end_tag, end_len) != 0)
      {      
         if(strstr(line, "Dimensions") != NULL)
         {
            pTok = &line[10];
            //extract values, one-by-one
            for(k = 0; k < num; k++)
            {
               j = ExtractString(pTok, tmp);
               j = ValidateExtraction(j, k, num, "GridAlgorithm::InitFromFile()");
               pTok += j;            
               m_pDims[k] = atoi(tmp);

               if(m_pDims[k] < 2)
               {
                  sprintf(tmp, "Invalid grid dimension (%d), defaulting to 2", m_pDims[k]);
                  LogError(ERR_FILE_IO, tmp);
                  m_pDims[k] = 2;
               }
            }
         }/* end for() */                  
         else if(strstr(line, "EvalsPerIter") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MiniSize);            
         }
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   fclose(pFile);

   //determine grid size
   m_GridSize = 1;
   for(i = 0; i < num; i++){ m_GridSize *= m_pDims[i];}
   if(m_GridSize < 1)
   {
      sprintf(tmp, "Invalid grid size: %d", m_GridSize);
      LogError(ERR_FILE_IO, tmp);
      ExitProgram(1);
   }
  
   //allocate mini grid
   NEW_PRINT("GridStruct", 1);
   m_pMini = new GridStruct;
   MEM_CHECK(m_pMini);
   if(m_MiniSize > m_GridSize){ m_MiniSize = m_GridSize;}
   m_pMini->nprm = num;
   NEW_PRINT("double", m_MiniSize);
   m_pMini->f  = new double[m_MiniSize];
   NEW_PRINT("double", num);
   m_pMini->dp = new double[num];
   NEW_PRINT("double *", m_MiniSize);
   m_pMini->p  = new double*[m_MiniSize];
   for(i = 0; i < m_MiniSize; i++)
   { 
      NEW_PRINT("double", num);
      m_pMini->p[i]  = new double[num];
   }
   MEM_CHECK(m_pMini->p[i-1]);

   /* Determine number of iterations required (each iteration 
   will evaluate a single mini-grid) */
   m_NumIters = m_GridSize / m_MiniSize;
   if((m_GridSize % m_MiniSize) != 0) m_NumIters++;

   //initialize dp (step size for each parameter)
   for(i = 0; i < num; i++)
   { 
      dim = (double)m_pDims[i];
      upr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUprBnd();
      lwr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLwrBnd();
      m_pMini->dp[i] = (upr - lwr)/(dim - 1.00);
      m_Lwr[i] = lwr;
   }

   //allocate and compute the 'rollover' count for each parameter
   NEW_PRINT("int", num);
   m_Rval = new int[num];
   MEM_CHECK(m_Rval);
   
   m_Rval[num-1] = m_pDims[num-1];
   for(i = num-2; i >= 0; i--) { m_Rval[i] = m_Rval[i+1]*m_pDims[i];}
} /* end InitFromFile() */

/******************************************************************************
GetGridVal()

Compute value of parameter (j) at grid location (i).
******************************************************************************/
double GridAlgorithm::GetGridVal(int i, int j)
{
   static int last = m_pMini->nprm - 1;
   int idx;
   if(j == last) idx = (i % m_Rval[j]);
   else          idx = (i % m_Rval[j])/m_Rval[j+1];
   return(m_Lwr[j] + m_pMini->dp[j]*idx);
}/* end GetGridVal() */

/******************************************************************************
GRID_Program()

Calibrate or optimize the model using the grid algorithm.
******************************************************************************/
void GRID_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("GridAlgorithm", 1);
   GridAlgorithm * Grid = new GridAlgorithm(model);
   MEM_CHECK(Grid);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { Grid->Calibrate(); }
   else { Grid->Optimize(); }

   delete Grid;
   delete model;
} /* end GRID_Program() */
