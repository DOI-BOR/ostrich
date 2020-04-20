#include "mpi_stub.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "PDDSAlgorithm.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

//for keeping track of the state of the algorithm
#define DDS_INIT_STATE   (0)
#define DDS_SEARCH_STATE (1)
#define DDS_DONE_STATE   (2)

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void PDDSAlgorithm::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   m_NumInit = 1;
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
max_int()

Max value of two integers.
******************************************************************************/
int PDDSAlgorithm::max_int(int a, int b) 
{
   if(a > b) return a;
   return b;
}/* end max_int() */

/**********************************************************************
ResetUserSeed()

**********************************************************************/
void PDDSAlgorithm::ResetUserSeed(int seed)
{
   int i;
   m_UserSeed = seed; 
   ResetRandomSeed(seed);
   int NRmax = 10*(m_MaxIter*m_num_dec);
   for(i = 0; i < NRmax; i++)
   {
      m_harvest[i] = random_number();
   }
}/* end ResetUserSeed() */

/******************************************************************************
PDDSAlgorithm CTOR
******************************************************************************/
PDDSAlgorithm::PDDSAlgorithm(ModelABC *pModel)
{
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char begin_token[DEF_STR_SZ];
   char end_token[DEF_STR_SZ];
	
   RegisterAlgPtr(this);

   m_pModel = pModel;
   m_pStats = NULL;
	
   //init. everything to reasonable defaults
   m_r_val						= 0.2;
   m_alpha                 = 0.5;
   m_beta                  = 0.5;
   m_UserSeed				   = GetRandomSeed();
   m_CurIter               = 0;
   m_MaxIter					= 100;
   m_UserSuppliedInit = false;
   m_DEBUG_dds = false;
   m_DEBUG_neigh_value = false;
   strcpy(m_use_opt, "standard");

   //Read Data From Algorithm Input file 
   IroncladString pFileName = GetInFileName();
   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("PDDSAlgorithm::CTOR", pFileName);
   }/* end if() */ 

  //accept multiple section headings
  if(CheckToken(inFile, "BeginParallelDDSAlg", pFileName) == true)
  {
     strcpy(begin_token, "BeginParallelDDSAlg");
     strcpy(end_token, "EndParallelDDSAlg");
  }
  else
  {
     rewind(inFile);
     if(CheckToken(inFile, "BeginParallelDDS", pFileName) == true)
     {
        strcpy(begin_token, "BeginParallelDDS");
        strcpy(end_token, "EndParallelDDS");
     }
     else
     {
        rewind(inFile);
        if(CheckToken(inFile, "BeginParaDDSAlg", pFileName) == true)
        {
           strcpy(begin_token, "BeginParaDDSAlg");
           strcpy(end_token, "EndParaDDSAlg");
        }
        else
        {
           rewind(inFile);
           if(CheckToken(inFile, "BeginParaDDS", pFileName) == true)
           {
              strcpy(begin_token, "BeginParaDDS");
              strcpy(end_token, "EndParaDDS");
           }
           else
           {
               rewind(inFile);
               if(CheckToken(inFile, "BeginDDSAlg", pFileName) == true)
               {
                  strcpy(begin_token, "BeginDDSAlg");
                  strcpy(end_token, "EndDDSAlg");
               }
               else
               {
                  rewind(inFile);
                  if(CheckToken(inFile, "BeginDDS", pFileName) == true)
                  {
                     strcpy(begin_token, "BeginDDS");
                     strcpy(end_token, "EndDDS");
                  }
                  else //default
                  {
                     strcpy(begin_token, "BeginDDS");
                     strcpy(end_token, "EndDDS");
                  }
               }
            }
         }
      }
   }
   rewind(inFile);

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
        else if (strstr(line, "UseOpt") != NULL)
         {
            sscanf(line, "%s %s", tmp, &(m_use_opt[0]));
         }
          else if (strstr(line, "AlphaValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_alpha);
         }
         else if (strstr(line, "BetaValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_beta);
         }
         else if(strstr(line, "MaxIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxIter);
         }
         else if(strstr(line, "UseInitialParamValues") != NULL)
         {
            m_UserSuppliedInit=true;
         }
         else if(strstr(line, "UseRandomParamValues") != NULL)
         {
            m_UserSuppliedInit=false;
         }
         else if(strstr(line, "EnableDebugging") != NULL)
         {
            m_DEBUG_dds=true;
            m_DEBUG_neigh_value=true;
         }

         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default PDDS algorithm setup.");
   }/* end else() */

   if ((m_r_val<0) || (m_r_val>1))
   {
      LogError(ERR_FILE_IO,"Bad Perturbation value specified for DDS Algorithm");
      ExitProgram(1);
   }

   if (m_MaxIter<1)
   {
      LogError(ERR_FILE_IO,"Maximum DDS Algorithm iterations must be >0");
      ExitProgram(1);
   }

   /* read in one or more initial parameter sets specified by user */
   rewind(inFile);
   int i, j, k;
   int num = 0;
   m_NumInit = 0;
   m_pInit = NULL;
   char * pTok;
   if(CheckToken(inFile, "BeginInitParams", pFileName) == true)
   {
      FindToken(inFile, "EndInitParams", pFileName);
      rewind(inFile);

      //allocate space for the parameter list
      num = m_pModel->GetParamGroupPtr()->GetNumParams();

      //count the number of entries
      FindToken(inFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      m_NumInit = 0;
      while(strstr(line, "EndInitParams") == NULL)
      {
         m_NumInit++;
         line = GetNxtDataLine(inFile, pFileName);
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
      rewind(inFile);
      FindToken(inFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      i = 0;
      while(strstr(line, "EndInitParams") == NULL)
      {
         pTok = line;
         //extract values, one-by-one, making any necessary conversions
         for(k = 0; k < num; k++)
         {
            j = ExtractString(pTok, tmp);
            j = ValidateExtraction(j, k, num, "PDDS::CTOR()");
            pTok += j;            
            m_pInit[i][k] = m_pModel->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
         }/* end for() */                  
         i++;
         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   }/* end if() */
   fclose(inFile);

   InitDdsDataMembers();

   IncCtorCount();
}/* end CTOR() */

/**********************************************************************
Destroy()
**********************************************************************/
void PDDSAlgorithm::Destroy()
{
   int i;

   delete m_pStats;

   for(i = 0; i < m_NumInit; i++)
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;

   DestroyDdsDataMembers();
   IncDtorCount();
}/* end DTOPR () */

/**********************************************************************
Optimize()

Coverted to C from file = dds.f90_MPI

Dynamically dimensioned Search (DDS) version 1.1 algorithm by Bryan Tolson
Ostrich version from Fortran version (original was coded in Matlab)
Coded in Nov 2005 by Bryan Tolson

DDS is an n-dimensional continuous global optimization algorithm.  It is coded as a 
minimizer but built into the code is a hidden transformation to make it capable of 
solving a maximization problem.  In a maximization problem, the algorithm minimizes 
the negative of the objective function F (-1*F).  Ostrich ensures it is always a
min problem.
       
REFERENCE FOR THIS ALGORITHM:  
Tolson, B. A., and C. A. Shoemaker (2007), Dynamically dimensioned search algorithm 
for computationally efficient watershed model calibration, Water Resour. Res., 43, 
W01413, doi:10.1029/2005WR004723.
**********************************************************************/
void PDDSAlgorithm::Optimize(void)
{
   // IMPORTANT CLASS DEFINITIONS
   // s_min and s_max : decision variable bounds
   // maxiter : the total allowable objective function evaluations
   // r_val : perturbation size parameter -> 0.2 by default
   
   // IMPORTANT MEMBER FUNCTIONS
   // obj_func : user specified objective function (from Ostrich model class).
   // random_number : uniform random # generator
   // neigh_value : 1-dimensional decision variable perturbation routine
   int eval, ini_fevals, num_recv;
   int jct;
   int j, k, dvn_count, dv;
   int jj = 0;
   long int ileft;
   double ranval;
   double Ftest,fvalue,Pn,new_value;
   int state = DDS_INIT_STATE;
   bool DEBUG;
   bool bWarmStart;
   bool bSynch = SynchReceives();
   FILE * pPnFile;

   //variables for MPI/parallel 
   int number_of_times_slave_worked;
   bool work_left;
   int nslaves, nxtsid, slaveindex;
   const int tag = 0;
   MPI_Status status;

   int signal;
   const int dowork=101;
   const int stopwork=102;
   
   bWarmStart = m_pModel->CheckWarmStart();

   // maximum number of slaves, increase if needed
   int maxslaves = m_nprocessors+1;
   double ** slave_working_on_x;
   //allocate storage for slaves.
   slave_working_on_x = new double * [maxslaves];
   for(j = 0; j < maxslaves; j++)
   {
      slave_working_on_x[j] = new double[m_num_dec+1];
   }/* end for() */

   // set debug level
   DEBUG = m_DEBUG_dds;

   bool bBanner = false;
   double a = 0.00; //the fraction of elapsed budget
	StatusStruct    pStatus;

   //pre-emption variables
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();
   int nSpecial = pGroup->GetNumSpecialParams();
   double Fbest, * Cbest;
   
   //write setup
   WriteSetup(m_pModel, "Parallel Dynamically Dimensioned Search Algorithm (PDDS)");
   //write banner
   WriteBanner(m_pModel, "trial    best fitness   ", " trials remaining");
	pStatus.maxIter = m_MaxIter;

   /*--------------------------------------------------------------------------------------
   DDS *initialization* procedure: code in this section is not part of the DDS algorithm as 
   presented in Figure 1 of Tolson and Shoemaker (2007).  It was used in paper, just not 
   included in Figure 1 pseudocode.
   --------------------------------------------------------------------------------------*/
   //master section of code
   m_master = 0;
   num_recv = 0;
   ileft = m_MaxIter;
   if (m_rank == m_master)
   {
      Cbest = new double[nSpecial];
      for(int iS = 0; iS < nSpecial; iS++)
      {
         Cbest[iS] = 0.00;
      }

      /****************************************************************************
   
           This is the MPI version of DDS
   
      ****************************************************************************/

      pPnFile = fopen("OstDDSPn.txt", "w");
      fprintf(pPnFile, "EVAL  Pn\n");
      fclose(pPnFile);

      // Random Numbers in array harvest generated in DDS_allocate.f90
      m_ign = 1; // global variable used as index for harvest in neigh_value.f90

      /*------------------------------------------------------------------------------------
      start the OUTER DDS ALGORITHM LOOP for remaining allowble function evaluations (ileft)
      ------------------------------------------------------------------------------------*/
      jct = 0;

      nslaves = m_nprocessors - 1;
      slaveindex = 0;
      nxtsid = 0;

      for(eval = 1; eval <= (ileft + nslaves); eval++)
      {      
         //update status and check for user-selected quit
         pStatus.curIter = eval;
         if(IsQuit() == true){ MPI_Abort(MPI_COMM_WORLD,0); break;}

         // receive new f from slave, evaluate

         /*-----------------------------------------------------------------------
         If expecting a message from slave reporting back with work
         this will happen only after eval is iterated over the first nslaves iterations
         -----------------------------------------------------------------------*/
         if (eval > nslaves)
         {
            if(bSynch == true)
            {
               slaveindex = nxtsid + 1;
               nxtsid = (nxtsid + 1) % nslaves;
            }
            else
            {
               slaveindex = MPI_ANY_SOURCE;
            }
            /* -------------------------------------------------------------------
            slave will have done this: call obj_func(num_dec,stest,fvalue)
            obtain a evaluated solution from slave
            -------------------------------------------------------------------*/            
            MPI_Recv(&(m_stest[m_num_dec]),1+nSpecial,MPI_DOUBLE,slaveindex,tag,MPI_COMM_WORLD,&status); 
            fvalue = m_stest[m_num_dec];
            num_recv++;

            //determine source
            slaveindex = status.MPI_SOURCE;

            /*
            FILE * pLog = fopen("OstMessages0.txt", "a");
            fprintf(pLog, "SlaveIndex = %02d\n", slaveindex);
            fclose(pLog);
            */

            if (DEBUG == true) 
            {
               printf("fvalue = %E\n", fvalue);
            }/* end if() */
            Ftest = m_to_max * fvalue; // to_max handles min (=1) and max (=-1) problems, 
            
            if(bBanner == true)
            {
               WriteInnerEval(WRITE_DDS, 0, '.');
               bBanner = false;
               jj=0;
            }
                
            // update current best solution
            if ((num_recv == 1) || (Ftest <= m_Fbest))
            {               
               m_Fbest = Ftest;
               for(int iS = 0; iS < nSpecial; iS++)
               {
                  Cbest[iS] = m_stest[m_num_dec+1+iS];
               }
               if (DEBUG == true) 
               {
                  printf("%4d\t%E\n", eval, m_Fbest);
               }
               for(j = 0; j < m_num_dec; j++)
               {
                 m_sbest[j] = slave_working_on_x[slaveindex][j];         
               }/* end for() */

               /* --------------------------------------------------------
               Update status file for long runs so that LAST best-results
               are saved in status.out and status.bin
               -------------------------------------------------------- */
               FILE * pOut = fopen("dds_status.out", "a");               
               fprintf(pOut, "%4d\t%E",  num_recv, m_to_max*m_Fbest);
               for(k = 0; k < m_num_dec; k++)
               {
                  fprintf(pOut, "\t%E", m_sbest[k]);
               }
               fprintf(pOut, "\n");
               fclose(pOut);

               //Ostrich formatted output
               if(num_recv > 1)
               {
                  WriteInnerEval(++jj, 0, '.');
                  WriteInnerEval(WRITE_ENDED, 0, '.');
               }
               bBanner = true;
               pGroup->WriteParams(m_sbest);
               WriteRecord(m_pModel, num_recv, m_Fbest, (double)(m_MaxIter-num_recv));
               pStatus.pct  = ((float)(100)*(float)(num_recv)/(float)(m_MaxIter));
               pStatus.numRuns = m_pModel->GetCounter();
               WriteStatus(&pStatus);
               m_pModel->SaveBest(slaveindex); //save the input and output files of the best configuration              
            } /* end if() */
            else if(num_recv >= m_MaxIter)
            {
               WriteInnerEval(++jj, 0, '.');
               WriteInnerEval(WRITE_ENDED, 0, '.');
               pGroup->WriteParams(m_sbest);
               WriteRecord(m_pModel, num_recv, m_Fbest, 0.00);
            }
            else
            {
               WriteInnerEval(++jj, 0, '.');
            }
         }/* end if(i > nslaves) */
         else
         {
            /* ------------------------------------------------------------
            nothing to be done here except:
               loop over slaves to assign initial work
            ------------------------------------------------------------ */
            slaveindex = slaveindex + 1;
         }/* end else() */

         //if work remains to be done
         if(eval <= ileft)
         {
            /* -------------------------------------------------------------
            PAWEL this is the place to accumulate search history
            accumulate DDS search history
            --------------------------------------------------------------*/
            if(state == DDS_INIT_STATE)
            {
               //first time through, set number of initialization steps
               if(eval == 1) 
               {
                  if(bWarmStart == true)
                  {
                     // at least one function evaluation per slave used to initialize DDS
                     ini_fevals = max_int(nslaves, 1);
                     if(DEBUG == true)
                     {
                        printf("Warm start detected.... ");
                     }
                  }
                  else if (m_UserSuppliedInit == true)
                  {
                     // at least one function evaluation per slave used to initialize DDS
                     ini_fevals = max_int(nslaves, m_NumInit);
                     if(DEBUG == true)
                     {
                        printf("Evaluating user supplied initial solution.... ");
                     }
                  }/* end if() */
                  // standard random initialization as described in Tolson & Shoemaker (2007).
                  else
                  {
                     // use best of at least 5 (or nslaves) uniform random solutions.
                     // Calculate # function evals to use for uniform random sampling to initialize DDS 
                     ini_fevals = max_int(5,(int)(0.005*((double)m_MaxIter)));
                     if(ini_fevals < nslaves)
                     {
                        ini_fevals = nslaves;
                     }
                     if(DEBUG == true)
                     {
                        printf("Sampling for initial DDS solution....   ");
                     }
                  }/* end else() */

                  if (DEBUG == true)
                  {
                     printf("ini_fevals, ileft = %d , %ld\n", ini_fevals, ileft);
                  }

                  //initial banner in output file
                  FILE * pOut = fopen("dds_status.out", "w");
                  fprintf(pOut, "STEP\tOBJ._FUNCTION");
                  for(k = 0; k < m_num_dec; k++)
                  {
                     fprintf(pOut, "\t%-12s", m_DVnames[k]);
                  }
                  fprintf(pOut, "\n");
                  fclose(pOut);
               }/* end if(eval == 1) */

               if (DEBUG == true) 
               {
                  printf("i , ini_fevals = %d , %d\n", eval, ini_fevals);
               }

               /* ----------------------------------------------------------------------------------------
               Sample an initial solution candidate (uniform random sampling) unless there is a user 
               supplied solution available.
               ---------------------------------------------------------------------------------------- */
               if (((m_UserSuppliedInit == false) && (bWarmStart == false)) || (eval > m_NumInit))
               {
                  for(j = 0; j < m_num_dec; j++)
                  {
                     ranval = random_number();
                     m_stest[j] = m_s_min[j] + ranval * (m_s_max[j] - m_s_min[j]);
                  } /* end for() */
               } 
               else if (bWarmStart == true)
               {
                  WarmStart();
                  for(j = 0; j < m_num_dec; j++)
                  {
                     m_stest[j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal(); 
                  }
               }
               else
               {
                  for(j = 0; j < m_num_dec; j++)
                  {
                     m_stest[j] = m_pInit[eval-1][j]; // solution to evaluate is given by user
                  }
               }/* end if() */

               //advance to search phase when initialization is done
               if(eval == ini_fevals)
               {
                  state = DDS_SEARCH_STATE;
                  if(DEBUG == true)
                  {
                     printf("Done DDS initialization.\nDDS is running...\nStep\tFbest\n");
                  }
               }
            }
            else //(state == DDS_SEARCH_STATE)
            {
               /* ------------------------------------------------------------
               generate new trial value      
               Determine Decision Variable (DV) selected for perturbation:
               --------------------------------------------------------------*/
               if (strcmp(m_use_opt, "no-rand-num") == 0)
               {
                  Pn = m_alpha;
               }
               else
               {
                  if(eval <= 2*nslaves)
                  {
                     Pn = 1.00; //have each slave perturb all parameters the first time through
                  }
                  else
                  {
                     Pn = 1.0 - log((double)(eval-2*nslaves)) / log((double)(ileft-2*nslaves)); // probability each DV selected
                  }
                  if (DEBUG == true)
                  {
                     printf("Pn = %e\n", Pn);
                  }
               }
               pPnFile = fopen("OstDDSPn.txt", "a");
               fprintf(pPnFile, "%04d  %E\n", eval, Pn);
               fclose(pPnFile);
           
               dvn_count = 0; // counter for how many DVs selected for perturbation

               // define stest initially as best current solution
               for(j = 0; j < m_num_dec; j++)
               {
                  m_stest[j] = m_sbest[j];
               }/* end for() */
                
               if (DEBUG == true)
               {
                  printf("stest\n");
                  for(j = 0; j < m_num_dec; j++)
                  {
                     printf("%E\n", m_stest[j]);
                  }/* end for() */
               }/* end if() */
   
               for(j = 0; j < m_num_dec; j++)
               {
                  // Using ngd instead of idg=idg+1 should help when ORDERED is removed
                  m_ngd = j + (eval-1)*m_num_dec + jct;
                  ranval = m_harvest[m_ngd];

                  if (DEBUG == true)
                  {
                     printf("eval , ranval, m_ngd = %d , %E, %d\n", eval, ranval, m_ngd);
                  }
          
                  // jth DV selected for perturbation
                  if (ranval < Pn) 
                  {
                     dvn_count = dvn_count + 1;
                     // call 1-D perturbation function to get new DV value (new_value):
                     new_value = neigh_value(m_sbest[j], m_s_min[j], m_s_max[j], m_r_val);
                     // note that r_val is the value of the DDS r perturbation size parameter (0.2 by default)
                     m_stest[j] = new_value; // change relevant DV value in stest
                  }/* end if() */
               } /* end for (dv perturbation) */

               if (DEBUG == true)
               {
                  printf("dvn_count=%d i=%d Pn=%E\n", dvn_count, eval, Pn);
               }/* end if () */

               // no DVs selected at random, so select ONE
               if (dvn_count == 0) 
               {
                  ranval = m_harvest[m_ngd + 1];
                  if (DEBUG == true) 
                  {
                     printf("eval,ranval,m_ngd=%d, %E, %d\n",eval,ranval,m_ngd+1);
                  }
                  jct = jct + 1;
                  dv = (int)(ceil((double)(m_num_dec*ranval)))-1; // 0-based index for one DV
                  // call 1-D perturbation function to get new DV value (new_value):
                  new_value = neigh_value(m_sbest[dv],m_s_min[dv],m_s_max[dv],m_r_val);
                  m_stest[dv] = new_value; // change relevant DV value in stest
               }/* end if(pick ONE dv to perturb) */
            }/* end if(in search state) */

            // Tell slave to evaluate obj function value (fvalue) for stest, for example see grie10.f:
            if (DEBUG == true)
            {            
               printf("num_dec = %d\n",m_num_dec);
               printf("stest\n");
               for(j = 0; j < m_num_dec; j++)
               {
                  printf("%E\n", m_stest[j]);
               }/* end for() */
            }

            // send work to slave
            signal = dowork;
            MPI_Send(&signal,1,MPI_INT,slaveindex,tag,MPI_COMM_WORLD); 
            m_stest[m_num_dec] = m_Fbest;
            for(int iS = 0; iS < nSpecial; iS++)
            {
               m_stest[m_num_dec+1+iS] = Cbest[iS];
            }
            MPI_Send(m_stest,m_num_dec+1+nSpecial,MPI_DOUBLE,slaveindex,tag,MPI_COMM_WORLD);  
            for(j = 0; j < m_num_dec; j++)
            {
               slave_working_on_x[slaveindex][j] = m_stest[j];
            }/* end for() */
         }/* end if(work remaining) */
         else // no work remains, so send termination signals
         {
            state = DDS_DONE_STATE;
            signal = stopwork;
            MPI_Send(&signal,1,MPI_INT,slaveindex,tag,MPI_COMM_WORLD);
         }/* end else() */
      } /* end for(OUTER DDS ALGORITHM LOOP) */
   }/* end if(master) */
   else //PAWEL slave section of code
   {
      if(IsQuit() == true){ MPI_Abort(MPI_COMM_WORLD, 0);}

      number_of_times_slave_worked = 0;
      work_left = true;

      while(work_left == true)
      {
         MPI_Recv(&signal,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status); 
         if(signal == stopwork)
         {
            work_left = false;
            printf("termination signal received by process %d\n",m_rank);
         }
         else if(signal == dowork)
         {
            number_of_times_slave_worked = number_of_times_slave_worked + 1;
            MPI_Recv(m_stest,m_num_dec+nSpecial+1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status); 
			   //let special parameters know about best config.
            Fbest = m_stest[m_num_dec];
            Cbest = &(m_stest[m_num_dec+1]);
			   pGroup->ConfigureSpecialParams(Fbest, Cbest);

            //run model and let master know about reviesed special parameters
            m_stest[m_num_dec] = obj_func(m_num_dec,m_stest);
  			   pGroup->GetSpecialConstraints(Cbest);
            MPI_Send(&(m_stest[m_num_dec]),1+nSpecial,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);

            if(number_of_times_slave_worked == 1)
            {
               /* --------------------------------------------
               enable special parameters now that best 
               is initialized.
               -------------------------------------------- */
               pGroup->EnableSpecialParams();
            }
         }
         else
         {
            printf("unknown signal, error\n");
         }
      }/* end while() */

      printf("slave %d handled %d tasks\n", m_rank, number_of_times_slave_worked);
   } /* end if(slave section of code) */

   // PAWEL master section of code --- only master does output
   if(m_rank == m_master)
   {
	   pGroup->WriteParams(m_sbest); 
      m_pModel->Execute();
      WriteOptimal(m_pModel, m_Fbest);
      m_pModel->SaveBest(m_rank); //save the input and output files of the best configuration

	   pStatus.pct	= 100.0;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      WriteAlgMetrics(this);

      if(DEBUG == true)
      {
         printf("DDS optimization is complete!\n");
      }

      delete [] Cbest;
   }/* end if() */   
   return;
}/*end Optimize() */  

/*=========================================================================================
neigh_value()

   Ported to C from neigh_value.f90

   Purpose is to generate a neighboring decision variable value for a single
   decision variable value being perturbed by the DDS optimization algorithm.
   New DV value respects the upper and lower DV bounds.
   Coded by Bryan Tolson, Nov 2005.

   I/O variable definitions:
      x_cur - current decision variable (DV) value
      x_min - min DV value
      x_max - max DV value
      r  - the neighborhood perturbation factor
      new_value - new DV variable value (within specified min and max)
=========================================================================================*/
double PDDSAlgorithm::neigh_value(double x_cur, double x_min, double x_max, double r)
{
   double new_value;
   double ranval;
   double zvalue,x_range;
   double Work3, Work2, Work1;
   bool DEBUG;
 
   DEBUG = m_DEBUG_neigh_value;

   if (strcmp(m_use_opt, "no-rand-num") == 0)
   {
      new_value = x_cur * m_beta;
   }
   else
   {
      if (strcmp(m_use_opt, "try-int-solution") == 0)
      {
         new_value = (double)((int)(x_cur));
      }
      else
      {
         x_range = x_max - x_min;
         /*----------------------------------------------------------------------------------------
         Generate a standard normal random variate (zvalue). Perturb current value with normal 
         random variable. Below returns a standard Gaussian random number based upon well-known 
         Marsagalia-Bray Algorithm
         ----------------------------------------------------------------------------------------*/
         Work3=2.0;
         while ((Work3 >= 1.0) || (Work3 == 0.0))
         {
            ranval = m_harvest[m_ign];
            Work1 = 2.0 * (double)(ranval) - 1.0;
            ranval = m_harvest[m_ign+1];
            Work2 = 2.0 * (double)(ranval) - 1.0;
            Work3 = Work1 * Work1 + Work2 * Work2;
            m_ign = m_ign + 2;
         }/* end while() */  
         
         Work3 = sqrt((-2.0 * log(Work3)) / Work3); // natural log
        
         // pick one of two deviates at random (don't worry about trying to use both): 
         ranval = m_harvest[m_ign];
         m_ign = m_ign + 1;

         if (ranval <  0.5)
         {
            zvalue = Work1 * Work3;
         }
         else
         {
            zvalue = Work2 * Work3;
         }   
         // ------------ done standard normal random variate generation -------------------

         // calculate new decision variable value:                                
         new_value = x_cur + zvalue*r*x_range;
      }/* end else() */
   }/* end else() */

   // check new value is within DV bounds.  If not, bounds are reflecting.
   if (new_value < x_min)
   {
      new_value = x_min + (x_min - new_value);
      if (new_value > x_max)
      {
         /*--------------------------------------------------------------------
         If reflection goes past x_max then value should be x_min since 
         without reflection the approach goes way past lower bound.  
         This keeps x close to lower bound when x_cur is close to lower bound
         Practically speaking, this should never happen with r values <0.3
         ---------------------------------------------------------------------*/
         new_value = x_min;
      }
   }
   else if (new_value > x_max)
   {
      new_value= x_max - (new_value - x_max);
      if (new_value < x_min)
      {
         //see reasoning above 
         new_value = x_max;
      }
   }

   if(m_DEBUG_neigh_value == true)
   {
      printf("neigh_value = %f, index = %d\n", new_value, m_ign);
   }
   return new_value;
}/* end neigh_value() */

/*=========================================================================================
InitDdsDataMembers()

   Configure the DDS member variables.
=========================================================================================*/
void PDDSAlgorithm::InitDdsDataMembers(void)
{
   int i;

   m_master = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &(m_rank));
   MPI_Comm_size(MPI_COMM_WORLD, &(m_nprocessors));

   // DDS user input values extracted from Ostrich internals
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();
   m_num_dec = pGroup->GetNumParams();

   //DDS internal counters?
   m_ngd = 0;
   m_ign = 0;   

   // DDS Output declarations
   m_to_max = +1;

   //allocate storage for names of design vars
   m_DVnames = new char *[m_num_dec];
   for(i = 0; i < m_num_dec; i++)
   {
      m_DVnames[i] = new char[100];
      strncpy(m_DVnames[i], pGroup->GetParamPtr(i)->GetName(), 100);
   }/* end for() */

   //allocate 1-based storage for various design var vectors
   m_s_min = new double[m_num_dec];
   m_s_max = new double[m_num_dec];
   m_sbest = new double[m_num_dec];
   int nSpecial = pGroup->GetNumSpecialParams();
   m_stest = new double[m_num_dec+nSpecial+1];
   for(i = 0; i < m_num_dec; i++)
   {
      m_s_min[i] = pGroup->GetParamPtr(i)->GetLwrBnd();
      m_s_max[i] = pGroup->GetParamPtr(i)->GetUprBnd();
      m_sbest[i] = pGroup->GetParamPtr(i)->GetEstVal();
      m_stest[i] = m_sbest[i];
   }/* end for() */

   int NRmax = 10*(m_MaxIter*m_num_dec);
   m_harvest = new double[NRmax];
   for(i = 0; i < NRmax; i++)
   {
      m_harvest[i] = random_number();
   }
}/* end InitDdsDataMembers() */

/*=========================================================================================
DestroyDdsDataMembers()

   Free up dynamically allocated portions of the DDS class.
=========================================================================================*/
void PDDSAlgorithm::DestroyDdsDataMembers(void)
{
   int i;
   //free up storage for names of design vars
   for(i = 0; i < m_num_dec; i++)
   {
      delete [] m_DVnames[i];
   }/* end for() */
   delete [] m_DVnames;

   //free up 0-based storage for various design var vectors
   delete [] m_harvest;
   delete [] m_s_min;
   delete [] m_s_max;
   delete [] m_stest;   
   delete [] m_sbest;   
}/* end DestroyDdsDataMembers() */

/**********************************************************************
random_number()

	Return uniformly distributed random number between 0 and 1
   Coded by James Craig
**********************************************************************/
double PDDSAlgorithm::random_number(void)
{
	return (double)(MyRand()) / (double)(MY_RAND_MAX);
}/* end random_number() */

/******************************************************************************
obj_func()

User-defined objective function. A pass-through to the model class.
******************************************************************************/
double PDDSAlgorithm::obj_func(int nopt, double * x_values)
{
   static double a = 0.00;
   double fvalue;

   //correct parameter values if options are enabled
   MakeParameterCorrections(x_values, m_sbest, nopt, a);

   //pass final design vars along to parameter group class via model
   m_pModel->GetParamGroupPtr()->WriteParams(x_values);

   //run the model and capture result   
   m_pModel->Execute();
   m_CurIter++;

   a += 1.00/(double)m_MaxIter;

   fvalue = m_pModel->GetObjFuncVal();

   //send result back to DDS algorithm
   return fvalue;
}/* end obj_func() */

/*************************************************************************************
MakeParameerCorrections()
*************************************************************************************/
void PDDSAlgorithm::MakeParameterCorrections(double * x, double * xb, int n, double a)
{
   double lwr, upr;
   ParameterGroup * pParamGroup;
	pParamGroup = m_pModel->GetParamGroupPtr(); 

   for(int k = 0; k < n; k++)
   {
      lwr=pParamGroup->GetParamPtr(k)->GetLwrBnd();
      upr=pParamGroup->GetParamPtr(k)->GetUprBnd();
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
Calibrate()

Solve the Least-Squares minimization problem using PDDS
******************************************************************************/
void PDDSAlgorithm::Calibrate(void)
{ 
	char fileName[DEF_STR_SZ];
	FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);

	Optimize();

  sprintf(fileName, "OstOutput%d.txt", 0);

  //compute statistics (variance and covariance)
  m_pStats->CalcStats();

  //write statistics of best parameter set to output file
  pFile = fopen(fileName, "a");   
  m_pStats->WriteStats(pFile);
  fclose(pFile);

  //write statistics of best parameter set to output file
  m_pStats->WriteStats(stdout);

} /* end Calibrate() */

/**********************************************************************
WriteMetrics()
**********************************************************************/
void PDDSAlgorithm::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Parallel Dynamically-Dimensioned Search Algorithm (PDDS)\n");
   fprintf(pFile, "Desired Convergence Val : N/A\n");
   fprintf(pFile, "Actual Convergence Val  : N/A\n");
   fprintf(pFile, "Max Generations         : %d\n", m_MaxIter);
   fprintf(pFile, "Actual Generations      : %d\n", m_MaxIter);
   fprintf(pFile, "Peterbation Value       : %lf\n", m_r_val);

   fprintf(pFile, "Debug Statements        : ");
   if(m_DEBUG_dds == true)
   {
      fprintf(pFile, "enabled\n");
   }
   else
   {
      fprintf(pFile, "disabled\n");
   }

   fprintf(pFile, "Initial Solution        : ");
   if(m_UserSuppliedInit == true)
   {
      fprintf(pFile, "User Supplied\n");
   }
   else
   {
      fprintf(pFile, "Randomly Generated\n");
   }

   fprintf(pFile, "Special User Option     : %s\n", m_use_opt);
   if(strcmp(m_use_opt, "no-rand-num") == 0)
   {
      fprintf(pFile, "Alpha Value             : %lf\n", m_alpha);
      fprintf(pFile, "Beta Value              : %lf\n", m_beta);
   }
   else
   {
      fprintf(pFile, "Alpha Value             : not used\n");
      fprintf(pFile, "Beta Value              : not used\n");
   }

   fprintf(pFile, "Number of Processors    : %d\n", m_nprocessors);

   m_pModel->WriteMetrics(pFile);
   fprintf(pFile, "Algorithm successfully converged on a solution, however more runs may be needed\n");
}/* end WriteMetrics() */

/******************************************************************************
PDDS_Program()
Calibrate the model using PDDS.
******************************************************************************/
void PDDS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("DDS", 1);
   PDDSAlgorithm * PDDS = new PDDSAlgorithm(model);

   MEM_CHECK(PDDS);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { PDDS->Calibrate(); }
   else { PDDS->Optimize(); }

   delete PDDS;

   delete model;
} /* end PDDS_Program() */
