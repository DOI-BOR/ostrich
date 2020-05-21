//DDSAU.cpp
#include "mpi_stub.h"
#include <string.h>

#include "DDSAU.h"
#include "DDSAlgorithm.h"
#include "PDDSAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/**********************************************************************
DDSAU Constructor
**********************************************************************/
DDSAU::DDSAU(ModelABC *pModel)
{
	FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char yesno[DEF_STR_SZ];
   char begin_token[DEF_STR_SZ];
   char end_token[DEF_STR_SZ];
   int itmp;
	
	RegisterAlgPtr(this);

	m_pModel = pModel;
   m_pStats = NULL;
	
	//init. everything to reasonable defaults
   m_r_val = 0.2;
   m_nsols = 25;
   m_nbhvr = 0;
   m_MinIter = 30;
	m_MaxIter = 70;
   m_bParallel = false;
   m_fmax = 1000.00;
   m_bRandomize = false;
   m_bReviseAU = false;

	//Read Data From Algorithm Input file 
	IroncladString pFileName = GetInFileName();
	inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("DDSAU::CTOR", pFileName);
   }/* end if() */ 

  //set section tags
  strcpy(begin_token, "Begin_DDSAU_Alg");
  strcpy(end_token, "End_DDSAU_Alg");

  if(CheckToken(inFile, begin_token, pFileName) == true)
  {
      FindToken(inFile, end_token, pFileName);
      rewind(inFile);

      FindToken(inFile, begin_token, pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, end_token) == NULL)
      {
         if (strstr(line, "PerturbationValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_r_val);
         }
         else if (strstr(line, "NumSearches") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_nsols);
         }
         if (strstr(line, "Threshold") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_fmax);
         }
         else if(strstr(line, "MinItersPerSearch") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MinIter);
         }
         else if(strstr(line, "MaxItersPerSearch") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxIter);
         }
         else if(strstr(line, "ParallelSearches") != NULL)
         {
            sscanf(line, "%s %s", tmp, yesno);
            MyStrLwr(yesno);
            if(strcmp(yesno, "yes") == 0)
            {
               m_bParallel = true;
            }
         }
         else if(strstr(line, "Randomize") != NULL)
         {
            sscanf(line, "%s %s", tmp, yesno);
            MyStrLwr(yesno);
            if(strcmp(yesno, "yes") == 0)
            {
               m_bRandomize = true;
            }
         }
         else if(strstr(line, "ReviseAU") != NULL)
         {
            sscanf(line, "%s %s", tmp, yesno);
            MyStrLwr(yesno);
            if(strcmp(yesno, "no") == 0)
            {
               m_bReviseAU = false;
            }
         }
         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default DDSAU algorithm setup.");
   }/* end else() */

	 if ((m_r_val<0) || (m_r_val>1)){
     LogError(ERR_FILE_IO,"Bad Perturbation value specified for DDSAU Algorithm");
		 ExitProgram(1);
	 }
	 if (m_MaxIter<1){
     LogError(ERR_FILE_IO,"Maximum DDSAU Algorithm iterations must be >0");
		 ExitProgram(1);
	 }
	 if (m_MinIter<1){
     LogError(ERR_FILE_IO,"Minimum DDSAU Algorithm iterations must be >0");
		 ExitProgram(1);
	 }
    if (m_MaxIter<m_MinIter)
    {
      itmp = m_MaxIter;
      m_MaxIter = m_MinIter;
      m_MinIter = itmp;
    }

   fclose(inFile);

   SetObjFuncThreshold(m_fmax);

  m_pBehavioral = NULL;
  m_fBehavioral = NULL;
 
  IncCtorCount();
}

/**********************************************************************
		Destroy
**********************************************************************/
void   DDSAU::Destroy()
{
   int i;

   delete m_pStats;

   if(m_pBehavioral != NULL)
   {
      for(i = 0; i < m_nsols; i++)
      {
         delete [] m_pBehavioral[i];
      }
   }
   delete [] m_pBehavioral;
   delete [] m_fBehavioral;

  IncDtorCount();
}

/**********************************************************************
Optimize
**********************************************************************/
void   DDSAU::Optimize() 
{
   if(m_bParallel == false)
   {
      OptimizeSerial();
   }
   else
   {
      OptimizeParallel();
   }
}/* end Optimize() */

/**********************************************************************
OptimizeSerial()
**********************************************************************/
void   DDSAU::OptimizeSerial() 
{
   FILE * pOut, * pStdOut;
   DDSAlgorithm * pDDS;
   char outFileName[DEF_STR_SZ];
   int i, j, p, budget, nfx;
   char * line;
   int max_line_size;
   char * tmpStr;
   char * bestLine;
   int nBehavioral;
   int whichBehavioral;
   double fbest = 0.00;
   double fx;
   int iter;
   bool bOutFileExists;
   FILE * pTestFile;

   //write setup
   pStdOut = fopen("OstOutputDDSAU.txt", "w");
   fprintf(pStdOut, "DDS for Approximation of Uncertainty (DDSAU)\n");   
   //write banner


   fprintf(pStdOut, "Iter  Run   obj.function  ");
   for(p = 0; p < m_pModel->GetParamGroupPtr()->GetNumParams(); p++)
   {
      fprintf(pStdOut, "%-12s  ", m_pModel->GetParamGroupPtr()->GetParamPtr(p)->GetName());
   }
   fprintf(pStdOut, "\n");

   //perform the requested number of searches  
   pDDS = new DDSAlgorithm(m_pModel);
   for(i = 0; i < m_nsols; i++)
   {
      SetTrialNumber(i); //used by PreserveModelOutput

      /* ----------------------------
      Configure DDS search
      ---------------------------- */
      pDDS->SetPerturbationValue(m_r_val);
      pDDS->ResetUserSeed(1+GetRandomSeed());  //pass along a new random seed
      pDDS->SetNoUserInit();  //disable user init
      if(m_MinIter == m_MaxIter)
      {
         budget = m_MaxIter;
      }
      else
      {
         budget = m_MinIter + (MyRand() % ((m_MaxIter - m_MinIter) + 1));
      }
      pDDS->SetBudget(budget);

      /* -----------------------------------------------------------------
      Run DDS search
      ----------------------------------------------------------------- */

      //are there results from previous run?
      sprintf(outFileName, "OstModel0_DDS%d.txt", i);
      pTestFile = fopen(outFileName, "r");
      if(pTestFile == NULL)
      {
         bOutFileExists = false; 
      }
      else
      {
         fclose(pTestFile);
         bOutFileExists = true;
      }

      //if not revising, delete previous files and then optimize
      if(m_bReviseAU == false)
      {
         remove(outFileName);
         sprintf(outFileName, "OstOutput0_DDS%d.txt", i);
         remove(outFileName);
         pDDS->Optimize();
      }
      //revising previous results, but none found
      else if(bOutFileExists == false)
      {
         sprintf(outFileName, "OstOutput0_DDS%d.txt", i);
         remove(outFileName);
         pDDS->Optimize();
      }
      //re-use existing output files, no model runs required
      else 
      {
         printf("Using previous results located in %s\n", outFileName);

         remove("OstModel0.txt");
         rename(outFileName, "OstModel0.txt");

         remove("OstOutput0.txt");
         sprintf(outFileName, "OstOutput0_DDS%d.txt", i);
         rename(outFileName, "OstOutput0.txt");
      }

      /* ----------------------------
      Post-process DDS search
      ---------------------------- */
      max_line_size = GetMaxLineSizeInFile((char *)"OstModel0.txt");
      // space for lines of OstModel0.txt
      line = new char[max_line_size+1];
      line[0] = NULLSTR;
      tmpStr = new char[max_line_size+1];
      tmpStr[0] = NULLSTR;
      bestLine = new char[max_line_size+1];
      bestLine[0] = NULLSTR;
      //space for behavioral solutions
      if (m_pBehavioral == NULL)
      {
         m_pBehavioral = new char *[m_nsols];
         m_fBehavioral = new double[m_nsols];
         for (j = 0; j < m_nsols; j++)
         {
            m_pBehavioral[j] = new char[max_line_size+1];
            strcpy(m_pBehavioral[j], "unknown");
            m_fBehavioral[j] = NEARLY_HUGE;
         }
      }/* end if() */
      pOut = fopen("OstModel0.txt", "r");
      nfx = 0;
      nBehavioral = 0;
      
      fgets(line, max_line_size, pOut); //skip header
      while(!feof(pOut))
      {
         fgets(line, max_line_size, pOut);
         //is it a behavioral soloution?
         sscanf(line, "%d %lf %s", &iter, &fx, tmpStr);
         //count number of behavioral samples
         if(fx <= m_fmax)
         {
            nBehavioral++; 
         }/* end if() */
         //track best solution
         if((nfx == 0) || (fx < fbest))
         {
            fbest = fx;
            strcpy(bestLine, line);
         }
         nfx++;
      }/* end while() */
      fclose(pOut);

      //use the best solution
      if((m_bRandomize == false) || (nBehavioral == 0))
      {
         strncpy(m_pBehavioral[i], bestLine, max_line_size);
         m_fBehavioral[i] = fbest;
      }
      //randomly select a behavioral solution
      else
      {
         whichBehavioral = 1 + (MyRand() % nBehavioral);

         pOut = fopen("OstModel0.txt", "r");
         nBehavioral = 0;
         fgets(line, max_line_size, pOut); //skip header
         while(!feof(pOut))
         {
            fgets(line, max_line_size, pOut);
            //is it a behavioral soloution?
            sscanf(line, "%d %lf %s", &iter, &fx, tmpStr);
            //count number of behavioral samples
            if(fx <= m_fmax)
            {
               nBehavioral++; 
            }/* end if() */
            if(nBehavioral == whichBehavioral)
            {
               strcpy(m_pBehavioral[i], line);
               m_fBehavioral[i] = fx;
               break;
            }/* end if() */
         }/* end while() */
         fclose(pOut);         
      }/* end else() */

      if(m_fBehavioral[i] <= m_fmax)
      {
         m_nbhvr++;
      }/* end if() */

      /* -----------------------------------------------------------
      Results for each search will be stored in output files named 
      OstModel[N]_DDS[M].txt and OstOutput[N]_DDS[M].txt, where [N] 
      is the processor number and [M] is the DDS search number.
      ----------------------------------------------------------- */
      sprintf(outFileName, "OstModel0_DDS%d.txt", i);
      rename("OstModel0.txt", outFileName);
      sprintf(outFileName, "OstOutput0_DDS%d.txt", i);
      rename("OstOutput0.txt", outFileName);

      fprintf(pStdOut, "%-4d  %s", i, m_pBehavioral[i]);

      ((Model *)m_pModel)->SetCounter(0);

      delete [] line;
      delete [] tmpStr;
      delete [] bestLine;
   }/* end for() */		

   WriteMetrics(pStdOut);

   fclose(pStdOut);
   pDDS->Destroy();  
}/* end OptimizeSerial() */


/**********************************************************************
OptimizeParallel()
**********************************************************************/
void   DDSAU::OptimizeParallel() 
{
   FILE * pOut, * pStdOut;
   PDDSAlgorithm * pDDS;
   char outFileName[DEF_STR_SZ];
   char tmpFileName[DEF_STR_SZ];
   int i, p, budget, nfx, proc;
   char * line;
   char * tmpStr;
   char * bestLine;
   int nBehavioral;
   int whichBehavioral;
   double fbest = 0.00;
   double fx;
   int iter;
   int rank, nprocs, nFound;
   int max_line_size, max_line_size_i, j;
   FILE * pTestFile;

   MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
   MPI_Comm_size(MPI_COMM_WORLD, &(nprocs));
   
   pStdOut = NULL;
   if(rank == 0)
   {
      //write setup
      pStdOut = fopen("OstOutputDDSAU.txt", "w");
      fprintf(pStdOut, "Parallel DDS for Approximation of Uncertainty (PDDSAU)\n");   
      //write banner
      fprintf(pStdOut, "Iter  Run   obj.function  ");
      for(p = 0; p < m_pModel->GetParamGroupPtr()->GetNumParams(); p++)
      {
         fprintf(pStdOut, "%-12s  ", m_pModel->GetParamGroupPtr()->GetParamPtr(p)->GetName());
      }
      fprintf(pStdOut, "\n");
   }

   //perform the requested number of searches  
   pDDS = new PDDSAlgorithm(m_pModel);
   for(i = 0; i < m_nsols; i++)
   {
      SetTrialNumber(i);

      /* ----------------------------
      Configure DDS search
      ---------------------------- */
      pDDS->SetPerturbationValue(m_r_val);
      pDDS->ResetUserSeed(nprocs+GetRandomSeed());  //pass along a new random seed
      pDDS->SetNoUserInit();  //disable user init
      if(m_MinIter == m_MaxIter)
      {
         budget = m_MaxIter;
      }
      else
      {
         budget = m_MinIter + (MyRand() % ((m_MaxIter - m_MinIter) + 1));
      }
      pDDS->SetBudget(budget);

      /* -----------------------------------------------------------------
      Run DDS search
      ----------------------------------------------------------------- */
      //are there results from previous run?
      nFound = 0;
      if(rank == 0)
      {
         for(proc = 0; proc < nprocs; proc++)
         {
            sprintf(outFileName, "OstModel%d_DDS%d.txt", proc, i);
            pTestFile = fopen(outFileName, "r");
            if(pTestFile != NULL)
            {
               fclose(pTestFile);
               nFound++;
            }
         }/* end for(each processor) */
      }/* end if(rank == 0) */
      MPI_Bcast(&nFound, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

      //if not revising, delete previous files and then optimize
      if(m_bReviseAU == false)
      {
         if(rank == 0)
         {
            for(proc = 0; proc < nprocs; proc++)
            {
               sprintf(outFileName, "OstModel%d_DDS%d.txt", proc, i);
               remove(outFileName);
               sprintf(outFileName, "OstOutput%d_DDS%d.txt", proc, i);
               remove(outFileName);
            }
         }
         pDDS->Optimize();
      }/* end if(not revising) */
      //revising previous results, but not all output files were found
      else if(nFound != nprocs)
      {
         if(rank == 0)
         {
            for(proc = 0; proc < nprocs; proc++)
            {
               sprintf(outFileName, "OstModel%d_DDS%d.txt", proc, i);
               remove(outFileName);
               sprintf(outFileName, "OstOutput%d_DDS%d.txt", proc, i);
               remove(outFileName);
            }
         }
         pDDS->Optimize();
      }/* end if(no previous output) */
      //re-use existing output files, no model runs required
      else 
      {
         if(rank == 0)
         {
            for(proc = 0; proc < nprocs; proc++)
            {
               sprintf(outFileName, "OstModel%d_DDS%d.txt", proc, i);
               sprintf(tmpFileName, "OstModel%d.txt", proc);
               printf("Using previous results located in %s\n", outFileName);
               remove(tmpFileName);
               rename(outFileName, tmpFileName);   

               sprintf(outFileName, "OstOutput%d_DDS%d.txt", proc, i);
               sprintf(tmpFileName, "OstOutput%d.txt", proc);
               remove(tmpFileName);
               rename(outFileName, tmpFileName);   
            }/* end for() */
         }/* end rank() */
      }/* end else() */

      /* ----------------------------
      Post-process DDS search
      ---------------------------- */
      if(rank == 0)
      {
         nfx = 0;
         nBehavioral = 0;

         max_line_size_i = max_line_size = 0;
         for (proc = 0; proc < nprocs; proc++)
         {
            sprintf(outFileName, "OstModel%d.txt", proc);
            max_line_size_i = GetMaxLineSizeInFile(outFileName);
            if (max_line_size_i > max_line_size)
            {
               max_line_size = max_line_size_i;
            }/* end if() */
         }/* end for() */
         // space for lines
         line = new char[max_line_size+1];
         line[0] = NULLSTR;
         tmpStr = new char[max_line_size+1];
         tmpStr[0] = NULLSTR;
         bestLine = new char[max_line_size+1];
         bestLine[0] = NULLSTR;
         //space for behavioral solutions
         if(m_pBehavioral == NULL)
         {
            m_pBehavioral = new char *[m_nsols];
            m_fBehavioral = new double[m_nsols];
            for (j = 0; j < m_nsols; j++)
            {
               m_pBehavioral[j] = new char[max_line_size+1];
               strcpy(m_pBehavioral[j], "unknown");
               m_fBehavioral[j] = NEARLY_HUGE;
            }
         }

         for(proc = 0; proc < nprocs; proc++)
         {
            sprintf(outFileName, "OstModel%d.txt", proc);
            pOut = fopen(outFileName, "r");
            fgets(line, max_line_size, pOut); //skip header
            while(!feof(pOut))
            {
               fgets(line, max_line_size, pOut);
               //is it a behavioral soloution?
               sscanf(line, "%d %lf %s", &iter, &fx, tmpStr);
               //count number of behavioral samples
               if(fx <= m_fmax)
               {
                  nBehavioral++; 
               }/* end if() */
               //track best solution
               if((nfx == 0) || (fx < fbest))
               {
                  fbest = fx;
                  strcpy(bestLine, line);
               }
               nfx++;
            }/* end while() */
            fclose(pOut);
         }/* end for(each proc) */

         //use the best solution
         if((m_bRandomize == false) || (nBehavioral == 0))
         {
            strcpy(m_pBehavioral[i], bestLine);
            m_fBehavioral[i] = fbest;
         }
         //randomly select a behavioral solution
         else
         {
            whichBehavioral = 1 + (MyRand() % nBehavioral);

            nBehavioral = 0;
            for(proc = 0; proc < nprocs; proc++)
            {
               sprintf(outFileName, "OstModel%d.txt", proc);
               pOut = fopen(outFileName, "r");
               
               fgets(line, max_line_size, pOut); //skip header
               while(!feof(pOut))
               {
                  fgets(line, max_line_size, pOut);
                  //is it a behavioral soloution?
                  sscanf(line, "%d %lf %s", &iter, &fx, tmpStr);
                  //count number of behavioral samples
                  if(fx <= m_fmax)
                  {
                     nBehavioral++; 
                  }/* end if() */
                  if(nBehavioral == whichBehavioral)
                  {
                     strcpy(m_pBehavioral[i], line);
                     m_fBehavioral[i] = fx;
                     break;
                  }/* end if() */
               }/* end while() */
               fclose(pOut);         

               if(nBehavioral == whichBehavioral)
               {
                  break;
               }
            }/* end for() */
         }/* end else() */

         if(m_fBehavioral[i] <= m_fmax)
         {
            m_nbhvr++;
         }/* end if() */

         /* -----------------------------------------------------------
         Results for each search will be stored in output files named 
         OstModel[N]_DDS[M].txt and OstOutput[N]_DDS[M].txt, where [N] 
         is the processor number and [M] is the DDS search number.
         ----------------------------------------------------------- */
         for(proc = 0; proc < nprocs; proc++)
         {
            sprintf(outFileName, "OstModel%d_DDS%d.txt", proc, i);
            sprintf(tmpFileName, "OstModel%d.txt", proc);
            rename(tmpFileName, outFileName);

            sprintf(outFileName, "OstOutput%d_DDS%d.txt", proc, i);
            sprintf(tmpFileName, "OstOutput%d.txt", proc);
            rename(tmpFileName, outFileName);            
         }/* end for() */
         fprintf(pStdOut, "%-4d  %s", i, m_pBehavioral[i]);

         delete [] line;
         delete [] bestLine;
         delete [] tmpStr;
      }/* end if(rank == 0) */

      ((Model *)m_pModel)->SetCounter(0);
   }/* end for() */		

   if(rank == 0)
   {
      WriteMetrics(pStdOut);
      fclose(pStdOut);
   }
   pDDS->Destroy();  
}/* end OptimizeParallel() */

/*************************************************************************************
MakeParameerCorrections()
*************************************************************************************/
void DDSAU::MakeParameterCorrections(double * x, double * xb, int n, double a)
{
   return;
}/* enad MakeParameterCorrections() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using the DDS
******************************************************************************/
void DDSAU::Calibrate(void)
{     
	Optimize();
} /* end Calibrate() */

/**********************************************************************
WriteMetrics()
**********************************************************************/
void DDSAU::WriteMetrics(FILE * pFile)
{
   int i, p;

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm                : DDS for Approximating Uncertainty (DDSAU)\n");
   fprintf(pFile, "Perturbation Value       : %0.2lf\n", m_r_val);
   fprintf(pFile, "Desired # of Samples     : %d\n", m_nsols);
   fprintf(pFile, "Actual # of Samples      : %d\n", m_nbhvr);
   fprintf(pFile, "Min DDS Evals per Sample : %d\n", m_MinIter);
   fprintf(pFile, "Max DDS Evals per Sample : %d\n", m_MaxIter);
   fprintf(pFile, "Behavioral Threshold     : %E\n", m_fmax);
   if(m_bRandomize == true)
   {
      fprintf(pFile, "Randomize samples?       : yes\n");
   }
   else
   {
      fprintf(pFile, "Randomize samples?       : no\n");
   }
   if(m_bReviseAU == true)
   {
      fprintf(pFile, "Revise Previous DDS AU?  : yes\n");
   }
   else
   {
      fprintf(pFile, "Revise Previous DDS AU?  : no\n");
   }

   fprintf(pFile, "\nList of Behavioral Solutions\n");
   fprintf(pFile, "Iter  Run   obj.function  ");
   for(p = 0; p < m_pModel->GetParamGroupPtr()->GetNumParams(); p++)
   {
      fprintf(pFile, "%12s  ", m_pModel->GetParamGroupPtr()->GetParamPtr(p)->GetName());
   }
   fprintf(pFile, "\n");
   for(i = 0; i < m_nsols; i++)
   {
      if(m_fBehavioral[i] <= m_fmax)
      {
         fprintf(pFile, "%-4d  %s", i, m_pBehavioral[i]);
      }
   }/* end for() */
}/* end WriteMetrics() */

/******************************************************************************
DDSAU_Program()
Calibrate the model using DDS.
**************************************s****************************************/
void DDSAU_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("DDSAU", 1);
   DDSAU * pDDSAU = new DDSAU(model);

   MEM_CHECK(pDDSAU);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { pDDSAU->Calibrate(); }
   else { pDDSAU->Optimize(); }

   delete pDDSAU;
   delete model;
} /* end DDSAU_Program() */

