/******************************************************************************
File     : Ostrich.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

Main program execution. Provides a text interface for the set of optimization
and gridding algorithms that make up the Ostrich program.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added PSO
03-24-04    lsm   added PSO-LevMar hybrid, added ISOFIT_BUILD option
11-07-05    lsm   added support for BGA, GRID, VSA and CSA programs
03-03-07    jrc   added DDS program
******************************************************************************/
#include "mpi_stub.h"
#include <stdio.h>
#include <string.h>

#include "OptMathClass.h"
#include "OptSearchClass.h"
#include "SAAlgorithm.h"
#include "VandSA.h"
#include "APPSO.h"
#include "SCEUA.h"
#include "ComboSA.h"
#include "BisectionAlgorithm.h"
#include "PowellAlgorithm.h"
#include "SamplingAlgorithm.h"
#include "SteepDescAlgorithm.h"
#include "FletchReevesAlgorithm.h"
#include "LevenbergAlgorithm.h"
#include "GridAlgorithm.h"
#include "StatsClass.h"
#include "GeneticAlgorithm.h"
#include "BinaryGA.h"
#include "ParticleSwarm.h"
#include "DDSAlgorithm.h"/*JRC*/
#include "PDDSAlgorithm.h"
#include "DiscreteDDSAlgorithm.h"
#include "StatUtility.h"
#include "QuadTree.h"
#include "GLUE.h"
#include "RejectionSampler.h"
#include "SMOOTH.h"
#include "PADDS.h"
#include "ParaPADDS.h"
#include "BEERS.h"
#include "DDSAU.h"

#include "Exception.h"
#include "Utility.h"
#include "IsoParse.h"

#ifdef ISOFIT_BUILD
int Ostrich(int argc, StringType argv[])
#else
int main(int argc, StringType argv[])
#endif
{
   ProgramType program;
   UnmoveableString pInFile = GetOstFileName();

   /*-------------------------------------------------------------------
   Windows uses a non-standard 3-digit exponent that messes up 
   applications that rely on standard 2-digit fixed formatting.
   --------------------------------------------------------------------*/
#ifdef WIN32
      _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

   //initialize time tracker
   GetElapsedTime();

   double tStart = GetElapsedTics();
   MPI_Init(&argc,&argv);
   double tEnd = GetElapsedTics();
   int id;
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   if(id == 0) printf("Starting up MPI required %lf seconds\n", (tEnd - tStart));

   SetOstExePath(argv[0]);
#ifndef ISOFIT_BUILD
   InitErrors();   
#endif

   //initialize input files (assume only one input file)
   strcpy(GetOstFileName(), "ostIn.txt");
   strcpy(GetExeDirName(), ".");
   InitDataLine(GetOstFileName());
   
   program = ReadProgramType();
   SetProgramType(program);

   //execute desired operation
   switch(program)
   {
      case(GA_PROGRAM) : 
      {
         GA_Program(argc, argv);
         break;
      }/* end case(GA_PROGRAM) */
      case(BGA_PROGRAM) : 
      {
         BGA_Program(argc, argv);
         break;
      }/* end case(BGA_PROGRAM) */
      case(GRID_PROGRAM) : 
      {
         GRID_Program(argc, argv);
         break;
      }/* end case(GRID_PROGRAM) */
      case(SA_PROGRAM) : 
      {
         SA_Program(argc, argv);
         break;
      }/* end case(SA_PROGRAM) */
      case(CSA_PROGRAM) : //combinatorial simulated annealing
      {
         CSA_Program(argc, argv);
         break;
      }/* end case(CSA_PROGRAM) */
      case(VSA_PROGRAM) : //vanderbilt-louie simulated annealing
      {
         VSA_Program(argc, argv);
         break;
      }/* end case(VCSA_PROGRAM) */
      case(PSO_PROGRAM) : 
      {
         PSO_Program(argc, argv);
         break;
      }/* end case(PSO_PROGRAM) */
      case(PSO_LEV_PROGRAM) : 
      {
         PSO_LEVMAR_Program(argc, argv);
         break;
      }/* end case(PSO_LEV_PROGRAM) */
      case(APPSO_PROGRAM) : 
      {
         APPSO_Program(argc, argv);
         break;
      }/* end case(APPSO_PROGRAM) */
      case(SCEUA_PROGRAM) : 
      {
         SCEUA_Program(argc, argv);
         break;
      }/* end case(SCEUA_PROGRAM) */
      case(LEV_PROGRAM) : 
      {
         LEV_Program(argc, argv);
         break;
      }/* end case(LEV_PROGRAM) */
      case(GMLMS_PROGRAM) : 
      {
         GMLMS_Program(argc, argv);
         break;
      }/* end case(GMLMS_PROGRAM) */
      case(POWL_PROGRAM) : 
      {
         PWL_Program(argc, argv);
         break;
      }/* end case(POWL_PROGRAM) */
      case(STEEP_PROGRAM) : 
      {
         STPDSC_Program(argc, argv);
         break;
      }/* end case(STEEP_PROGRAM) */
      case(FLRV_PROGRAM) : 
      {
         FLRV_Program(argc, argv);
         break;
      }/* end case(FLRV_PROGRAM) */
      case(BIS_PROGRAM) : 
      {
         BIS_Program(argc, argv);
         break;
      }/* end case(BIS_PROGRAM) */
      case(SMP_PROGRAM) : 
      {
         SMP_Program(argc, argv);
         break;
      }/* end case(SMP_PROGRAM) */
      case(STATS_PROGRAM) : 
      {
         STATS_Program(argc, argv);
         break;
      }/* end case(STATS_PROGRAM) */
      case(JACOBIAN_PROGRAM) : 
      {
         Jacobian_Program(argc, argv);
         break;
      }/* end case(JACOBIAN_PROGRAM) */
      case(HESSIAN_PROGRAM) : 
      {
         Hessian_Program(argc, argv);
         break;
      }/* end case(HESSIAN_PROGRAM) */
      case(GRADIENT_PROGRAM) : 
      {
         Gradient_Program(argc, argv);
         break;
      }/* end case(GRADIENT_PROGRAM) */
      case(EVAL_PROGRAM) : 
      {
         EVAL_Program(argc, argv);
         break;
      }/* end case(EVAL_PROGRAM) */
      case(UTIL_PROGRAM) : 
      {
		 ConvertToASCII();
         //STATS_TestFdist();
         //STATS_TestStudentDist();
         //STATS_TestStdNormDist();
         break;
      }/* end case(UTIL_PROGRAM) */
      case(DDS_PROGRAM) :/*JRC*/
      {
         DDS_Program(argc,argv);
			break;
      }
      case(DDSAU_PROGRAM) :
      {
         DDSAU_Program(argc,argv);
			break;
      }
      case(PDDS_PROGRAM) :
      {
         PDDS_Program(argc,argv);
			break;
      }
      case(DDDS_PROGRAM) :
      {
         DiscreteDDS_Program(argc,argv);
			break;
      }
      case(GLUE_PROGRAM) :
      {
         GLUE_Program(argc,argv);
			break;
      }
      case(RJSMP_PROGRAM) :
      {
         RJSMP_Program(argc,argv);
			break;
      }
      case(METRO_PROGRAM) :
      {
         METRO_Program(argc,argv); //Metropolis MCMC
			break;
      }
      case(SMOOTH_PROGRAM) : 
      {
         SMOOTH_Program(argc, argv);
         break;
      }/* end case(SMOOTH_PROGRAM) */
      case(PADDS_PROGRAM) : 
      {
         PADDS_Program(argc, argv);
         break;
      }/* end case(PADDS_PROGRAM) */
      case(PARA_PADDS_PROGRAM) : 
      {
         PARA_PADDS_Program(argc, argv);
         break;
      }/* end case(PARA_PADDS_PROGRAM) */
      case(BEERS_PROGRAM) : 
      {
         BEERS_Program(argc, argv);
         break;
      }/* end case(BEERS_PROGRAM) */
      case(QUIT_PROGRAM) : 
      default:
      {            
         break;
      }/* end case(QUIT_PROGRAM) */      
   }/* end switch() */   

   tStart = GetElapsedTics();
   MPI_Finalize();
   tEnd = GetElapsedTics();
   if(id == 0) printf("Shutting down MPI required %lf seconds\n", (tEnd - tStart));

   #ifndef ISOFIT_BUILD
      ExitProgram(0);
   #endif

   MyRandCleanup();

   return (0);
} /* end main() */
