/******************************************************************************
File     : ParaPADDS.cpp
Author   : L. Shawn Matott
Copyright: 2015, L. Shawn Matott

The PADDS (Pareto Archived Dynamically Dimensioned Search) algorithm is a 
multi-objective version of the DDS algorithm. It has been ported from a C++ 
implementation of Mohammadamin Jahanpour (mjahanpo@uwaterloo.ca).

Version History
05-07-15    lsm   created file
******************************************************************************/

/*------------------------------------------------------------------- 
The copyright, disclaimer, and citation given below affects the 
following list of subroutines:

   bool covers(double* cub, double * regLow);
   bool partCovers(double* cub, double * regUp);
   int containsBoundary(double* cub, double * regLow, int split);
   double getMeasure(double * regLow, double * regUp);
   int isPile(double* cub, double * regLow, double * regUp);
   double getMedian(double * bounds, int size);
   double computeTrellis(double * regLow, double * regUp, double * trellis);
   void stream(double * regionLow, double * regionUp, double ** points, int npoints, int split, double cover);

Copyright (c) 2006 Nicola Beume  This program is free software: you 
can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, 
either version 2 of the License, or (at your option) any later version. 
This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License for more details. You should have received a copy of the 
GNU General Public License along with this program.

--------------------------------------------------------------------- 

This program calculates the dominated hypervolume or S-metric of a set 
of d-dimensional points (d>=3). Please refer to the following publication 
for a description of the algorithm: 
   Nicola Beume and Guenter Rudolph. Faster S-Metric Calculation by 
   Considering Dominated Hypervolume as Klee's Measure Problem. 
   In: B. Kovalerchuk (ed.): Proceedings of the Second IASTED Conference 
   on Computational Intelligence (CI 2006), pp. 231-236. 
   ACTA Press: Anaheim, 2006. 

   Extended version published as: Technical Report of the Collaborative 
   Research Centre 531 'Computational Intelligence', CI-216/06, 
   ISSN 1433-3325. University of Dortmund, July 2006. 
-------------------------------------------------------------------*/

#define PARA_PADDS_DEBUG (0)

#include <mpi.h>
#include <string.h>

#include "ParaPADDS.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ObjectiveFunction.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
ParaPADDS::ParaPADDS(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pNonDom = NULL;
   m_pDom = NULL;
   m_NumNonDom = 0;
   m_NumDom = 0;
   m_maxiter = 0;
   m_CurIter = 0;
   m_num_dec = 0;
   m_num_objs = 0;
   m_Select_metric = 3; //make "exact" the default
   m_fraction1 = 0.0;
   m_dominance_flag = 0;
   m_seed = 0;
   m_dim = 0;
   m_dimension = 0;
   m_dSqrtDataNumber = 0;
   m_volume = 0;
   m_stest_flat = NULL;
   m_InitTaskMethod = INIT_TASK_URANDOM;
   m_HamedOffset = 0;

   #if(PARA_PADDS_DEBUG == 1)
      printf("Created ParaPADDS instance\n");
      fflush(stdout);
   #endif
   
   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by PADDS and it's member variables.
******************************************************************************/
void ParaPADDS::Destroy(void)
{
   ArchiveStruct * pCur, * pDel;

   for(pCur = m_pNonDom; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   for(pCur = m_pDom; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   delete [] m_stest_flat;

   #if(PARA_PADDS_DEBUG == 1)
      printf("Destroyed ParaPADDS instance\n");
      fflush(stdout);
   #endif

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void ParaPADDS::WriteMetrics(FILE * pFile) 
{
   char select_str[DEF_STR_SZ];
   //0: Random
   //1: Crowding distance
   //2: Hypervolume Contribution (ESTIMATE)
   //3: Hypervolume Contribution (EXACT)
   if(m_Select_metric == 0) strcpy(select_str, "random");
   else if(m_Select_metric == 1) strcpy(select_str, "crowding distance");
   else if(m_Select_metric == 2) strcpy(select_str, "estimated hypervolume contribution");
   else if(m_Select_metric == 3) strcpy(select_str, "exact hypervolume contribution");
   else strcpy(select_str, "unknown");

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : PADDS - Pareto Archived Dynamically Dimensioned Search\n");
   fprintf(pFile, "Max Iterations          : %d\n", m_maxiter);
   fprintf(pFile, "Actual Iterations       : %d\n", m_CurIter);
   fprintf(pFile, "Num Decision Variables  : %d\n", m_num_dec);
   fprintf(pFile, "Num Objectives          : %d\n", m_num_objs);
   fprintf(pFile, "Random Seed             : %d\n", m_seed);
   fprintf(pFile, "Perturbation Value      : %lf\n", m_fraction1);
   fprintf(pFile, "Non-Dominated Solutions : %d\n", m_NumNonDom);  
   fprintf(pFile, "Dominated Solutions     : %d\n", m_NumDom);     
   fprintf(pFile, "Selection Metric        : %s\n", select_str);
   if(m_InitTaskMethod == INIT_TASK_URANDOM)
   {
      fprintf(pFile, "Initialization Method   : Uniform Random\n");
   }
   else if (m_InitTaskMethod == INIT_TASK_HAMED)
   {
      fprintf(pFile, "Initialization Method   : Hamed (Offset = %d)\n", m_HamedOffset);
   }
   else
   {
      fprintf(pFile, "Initialization Method   : Uniform Random\n");
   }

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void ParaPADDS::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];
   char begin_token[DEF_STR_SZ];
   char end_token[DEF_STR_SZ];
   
   m_maxiter = 50;
   m_fraction1 = 0.2;
   m_Select_metric = 3; //make "exact" the default

   //read in PADDS configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open PADDS config. file. Using Defaults");      
      return;
   }/* end if() */   

  //accept multiple section headings
  if(CheckToken(pFile, "BeginParallelPADDSAlg", pFileName) == true)
  {
     strcpy(begin_token, "BeginParallelPADDSAlg");
     strcpy(end_token, "EndParallelPADDSAlg");
  }
  else
  {
     rewind(pFile);
     if(CheckToken(pFile, "BeginParallelPADDS", pFileName) == true)
     {
        strcpy(begin_token, "BeginParallelPADDS");
        strcpy(end_token, "EndParallelPADDS");
     }
     else
     {
        rewind(pFile);
        if(CheckToken(pFile, "BeginParaPADDSAlg", pFileName) == true)
        {
           strcpy(begin_token, "BeginParaPADDSAlg");
           strcpy(end_token, "EndParaPADDSAlg");
        }
        else
        {
           rewind(pFile);
           if(CheckToken(pFile, "BeginParaPADDS", pFileName) == true)
           {
              strcpy(begin_token, "BeginParaPADDS");
              strcpy(end_token, "EndParaPADDS");
           }
           else
           {
               rewind(pFile);
               if(CheckToken(pFile, "BeginPADDSAlg", pFileName) == true)
               {
                  strcpy(begin_token, "BeginPADDSAlg");
                  strcpy(end_token, "EndPADDSAlg");
               }
               else
               {
                  rewind(pFile);
                  if(CheckToken(pFile, "BeginPADDS", pFileName) == true)
                  {
                     strcpy(begin_token, "BeginPADDS");
                     strcpy(end_token, "EndPADDS");
                  }
                  else //default
                  {
                     strcpy(begin_token, "BeginPADDS");
                     strcpy(end_token, "EndPADDS");
                  }
               }
            }
         }
      }
   }
   rewind(pFile);

   //make sure correct tokens are present
   if(CheckToken(pFile, begin_token, pFileName) == true)
   {
      FindToken(pFile, end_token, pFileName);
      rewind(pFile);

      FindToken(pFile, begin_token, pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, end_token) == NULL)
      {         
         if (strstr(line, "PerturbationValue") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_fraction1);
         }
         else if(strstr(line, "MaxIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_maxiter);
         }
         else if(strstr(line, "HamedOffset") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_HamedOffset);
         }
         else if(strstr(line, "SelectionMetric") != NULL)
         {
            sscanf(line, "SelectionMetric %s", tmp);

            //0: Random
            //1: Crowding distance
            //2: Hypervolume Contribution (ESTIMATE)
            //3: Hypervolume Contribution (EXACT)
            MyStrLwr(tmp);
            if(strcmp(tmp, "random") == 0) m_Select_metric = 0;
            else if(strcmp(tmp, "crowdingdistance") == 0) m_Select_metric = 1; 
            else if(strcmp(tmp, "estimatedhypervolumecontribution") == 0) m_Select_metric = 2; 
            else if(strcmp(tmp, "exacthypervolumecontribution") == 0) m_Select_metric = 3; 
            else m_Select_metric = 3; //make "exact" the default
         }
         else if(strstr(line, "InitializationTaskMethod") != NULL)
         {
            sscanf(line, "InitializationTaskMethod %s", tmp);
            // 0: Uniform random
            // 1: Hamed's Method
            MyStrLwr(tmp);
            if(strcmp(tmp, "urandom") == 0) m_InitTaskMethod = INIT_TASK_URANDOM;
            else if(strcmp(tmp, "hamed") == 0) m_InitTaskMethod = INIT_TASK_HAMED; 
            else m_InitTaskMethod = INIT_TASK_URANDOM; //make "urandom" the default
         }
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   fclose(pFile);

   #if(PARA_PADDS_DEBUG == 1)
      printf("Initialized ParaPADDS from file\n");
      printf("m_Select_metric = %d\n", m_Select_metric);
      printf("m_maxiter = %d\n", m_maxiter);
      printf("m_fraction1 = %lf\n", m_fraction1);
      fflush(stdout);
   #endif

} /* end InitFromFile() */

/******************************************************************************
InitUsingHamedMethod()

Initialize the solutions using Hamed's carefully constructed empirical method.
******************************************************************************/
void ParaPADDS::InitUsingHamedMethod(int n, double * stest, double * pMin, double * pMax, int offset, int sindex, int nsubs)
{
   // first pass - use urandom
   for (int j = 0; j < n; j++)
   {
      stest[j] = pMin[j] + (pMax[j] - pMin[j])*UniformRandom();
   }/* end for() */

   // second pass - update selected entries using Hamed's empirical method
   for (int j = 0; j < n/3; j++)
   {
      // The if structure distributes the initial conditions across the range of the "offset" term
      if (((double)sindex/(double)nsubs) <= 0.1) 
      {                           
         stest[j] = pMin[j];
         if ( (((double)sindex/(double)nsubs) >= 0.0) && (((double)sindex/(double)nsubs) <= 0.01) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.01) && (((double)sindex/(double)nsubs) <= 0.02)) {
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.02) && (((double)sindex/(double)nsubs) <= 0.03)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.03) && (((double)sindex/(double)nsubs) <= 0.04)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.04) && (((double)sindex/(double)nsubs) <= 0.05)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.05) && (((double)sindex/(double)nsubs) <= 0.06)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.06) && (((double)sindex/(double)nsubs) <= 0.07)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.07) && (((double)sindex/(double)nsubs) <= 0.08)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.08) && (((double)sindex/(double)nsubs) <= 0.09)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.09) && (((double)sindex/(double)nsubs) <= 0.098)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.098) && (((double)sindex/(double)nsubs) <= 0.10)){
            stest[j+offset] = pMax[j+offset]; }
     } 
     else if ( (((double)sindex/(double)nsubs) > 0.1) && (((double)sindex/(double)nsubs) <= 0.2) ) 
     {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.1;
         if ( (((double)sindex/(double)nsubs) > 0.10) && (((double)sindex/(double)nsubs) <= 0.11) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.11) && (((double)sindex/(double)nsubs) <= 0.12)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.12) && (((double)sindex/(double)nsubs) <= 0.13)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.13) && (((double)sindex/(double)nsubs) <= 0.14)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.14) && (((double)sindex/(double)nsubs) <= 0.15)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.15) && (((double)sindex/(double)nsubs) <= 0.16)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.16) && (((double)sindex/(double)nsubs) <= 0.17)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.17) && (((double)sindex/(double)nsubs) <= 0.18)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.18) && (((double)sindex/(double)nsubs) <= 0.19)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.19) && (((double)sindex/(double)nsubs) <= 0.198)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.198) && (((double)sindex/(double)nsubs) <= 0.20)){
            stest[j+offset] = pMax[j+offset]; }       
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.2) && (((double)sindex/(double)nsubs) <= 0.3) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.2;
         if ( (((double)sindex/(double)nsubs) > 0.20) && (((double)sindex/(double)nsubs) <= 0.21) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.21) && (((double)sindex/(double)nsubs) <= 0.22)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.22) && (((double)sindex/(double)nsubs) <= 0.23)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.23) && (((double)sindex/(double)nsubs) <= 0.24)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.24) && (((double)sindex/(double)nsubs) <= 0.25)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.25) && (((double)sindex/(double)nsubs) <= 0.26)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.26) && (((double)sindex/(double)nsubs) <= 0.27)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.27) && (((double)sindex/(double)nsubs) <= 0.28)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.28) && (((double)sindex/(double)nsubs) <= 0.29)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.29) && (((double)sindex/(double)nsubs) <= 0.298)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.298) && (((double)sindex/(double)nsubs) <= 0.30)){
            stest[j+offset] = pMax[j+offset]; }
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.3) && (((double)sindex/(double)nsubs) <= 0.4) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.3; 
         if ( (((double)sindex/(double)nsubs) > 0.30) && (((double)sindex/(double)nsubs) <= 0.31) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.31) && (((double)sindex/(double)nsubs) <= 0.32)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.32) && (((double)sindex/(double)nsubs) <= 0.33)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.33) && (((double)sindex/(double)nsubs) <= 0.34)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.34) && (((double)sindex/(double)nsubs) <= 0.35)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.35) && (((double)sindex/(double)nsubs) <= 0.36)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.36) && (((double)sindex/(double)nsubs) <= 0.37)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.37) && (((double)sindex/(double)nsubs) <= 0.38)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.38) && (((double)sindex/(double)nsubs) <= 0.39)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.39) && (((double)sindex/(double)nsubs) <= 0.398)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.398) && (((double)sindex/(double)nsubs) <= 0.40)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.4) && (((double)sindex/(double)nsubs) <= 0.5) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.4;
         if ( (((double)sindex/(double)nsubs) > 0.40) && (((double)sindex/(double)nsubs) <= 0.41) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.41) && (((double)sindex/(double)nsubs) <= 0.42)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.42) && (((double)sindex/(double)nsubs) <= 0.43)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.43) && (((double)sindex/(double)nsubs) <= 0.44)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.44) && (((double)sindex/(double)nsubs) <= 0.45)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.45) && (((double)sindex/(double)nsubs) <= 0.46)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.46) && (((double)sindex/(double)nsubs) <= 0.47)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.47) && (((double)sindex/(double)nsubs) <= 0.48)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.48) && (((double)sindex/(double)nsubs) <= 0.49)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.49) && (((double)sindex/(double)nsubs) <= 0.498)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.498) && (((double)sindex/(double)nsubs) <= 0.50)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.5) && (((double)sindex/(double)nsubs) <= 0.6) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.5;
         if ( (((double)sindex/(double)nsubs) > 0.50) && (((double)sindex/(double)nsubs) <= 0.51) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.51) && (((double)sindex/(double)nsubs) <= 0.52)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.52) && (((double)sindex/(double)nsubs) <= 0.53)){
             stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.53) && (((double)sindex/(double)nsubs) <= 0.54)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.54) && (((double)sindex/(double)nsubs) <= 0.55)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.55) && (((double)sindex/(double)nsubs) <= 0.56)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.56) && (((double)sindex/(double)nsubs) <= 0.57)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.57) && (((double)sindex/(double)nsubs) <= 0.58)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.58) && (((double)sindex/(double)nsubs) <= 0.59)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.59) && (((double)sindex/(double)nsubs) <= 0.598)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.598) && (((double)sindex/(double)nsubs) <= 0.60)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.6) && (((double)sindex/(double)nsubs) <= 0.7) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.6;
         if ( (((double)sindex/(double)nsubs) > 0.60) && (((double)sindex/(double)nsubs) <= 0.61) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.61) && (((double)sindex/(double)nsubs) <= 0.62)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.62) && (((double)sindex/(double)nsubs) <= 0.63)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.63) && (((double)sindex/(double)nsubs) <= 0.64)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.64) && (((double)sindex/(double)nsubs) <= 0.65)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.65) && (((double)sindex/(double)nsubs) <= 0.66)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.66) && (((double)sindex/(double)nsubs) <= 0.67)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.67) && (((double)sindex/(double)nsubs) <= 0.68)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.68) && (((double)sindex/(double)nsubs) <= 0.69)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.69) && (((double)sindex/(double)nsubs) <= 0.698)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.698) && (((double)sindex/(double)nsubs) <= 0.70)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.7) && (((double)sindex/(double)nsubs) <= 0.8) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.7;
         if ( (((double)sindex/(double)nsubs) > 0.70) && (((double)sindex/(double)nsubs) <= 0.71) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.71) && (((double)sindex/(double)nsubs) <= 0.72)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.72) && (((double)sindex/(double)nsubs) <= 0.73)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.73) && (((double)sindex/(double)nsubs) <= 0.74)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.74) && (((double)sindex/(double)nsubs) <= 0.75)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.75) && (((double)sindex/(double)nsubs) <= 0.76)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.76) && (((double)sindex/(double)nsubs) <= 0.77)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.77) && (((double)sindex/(double)nsubs) <= 0.78)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.78) && (((double)sindex/(double)nsubs) <= 0.79)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.79) && (((double)sindex/(double)nsubs) <= 0.798)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.798) && (((double)sindex/(double)nsubs) <= 0.80)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.8) && (((double)sindex/(double)nsubs) <= 0.9) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.8;
         if ( (((double)sindex/(double)nsubs) > 0.80) && (((double)sindex/(double)nsubs) <= 0.81) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.81) && (((double)sindex/(double)nsubs) <= 0.82)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.82) && (((double)sindex/(double)nsubs) <= 0.83)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.83) && (((double)sindex/(double)nsubs) <= 0.84)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.84) && (((double)sindex/(double)nsubs) <= 0.85)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.85) && (((double)sindex/(double)nsubs) <= 0.86)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.86) && (((double)sindex/(double)nsubs) <= 0.87)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.87) && (((double)sindex/(double)nsubs) <= 0.88)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.88) && (((double)sindex/(double)nsubs) <= 0.89)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.89) && (((double)sindex/(double)nsubs) <= 0.898)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.898) && (((double)sindex/(double)nsubs) <= 0.90)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.9) && (((double)sindex/(double)nsubs) <= 0.95) ) 
      {
         stest[j] = pMin[j] + (pMax[j] - pMin[j]) * 0.9;
         if ( (((double)sindex/(double)nsubs) > 0.90) && (((double)sindex/(double)nsubs) <= 0.91) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.91) && (((double)sindex/(double)nsubs) <= 0.92)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.25; 
         }else if((((double)sindex/(double)nsubs) > 0.92) && (((double)sindex/(double)nsubs) <= 0.93)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.50; 
         }else if((((double)sindex/(double)nsubs) > 0.93) && (((double)sindex/(double)nsubs) <= 0.94)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.75; 
         }else if((((double)sindex/(double)nsubs) > 0.94) && (((double)sindex/(double)nsubs) <= 0.95)){
            stest[j+offset] = pMax[j+offset]; }   
      } 
      else 
      {
         stest[j] = pMax[j];
         if ( (((double)sindex/(double)nsubs) > 0.95) && (((double)sindex/(double)nsubs) <= 0.96) ) {
            stest[j+offset] = pMin[j+offset]; 
         }else if((((double)sindex/(double)nsubs) > 0.96) && (((double)sindex/(double)nsubs) <= 0.97)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.25; 
         }else if((((double)sindex/(double)nsubs) > 0.97) && (((double)sindex/(double)nsubs) <= 0.98)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.50; 
         }else if((((double)sindex/(double)nsubs) > 0.98) && (((double)sindex/(double)nsubs) <= 0.99)){
            stest[j+offset] = pMin[j+offset] + (pMax[j+offset] - pMin[j+offset]) * 0.75; 
         }else if((((double)sindex/(double)nsubs) > 0.99) && (((double)sindex/(double)nsubs) <= 1)){
            stest[j+offset] = pMax[j+offset]; }
      }       
   } /* end for() */
} /* end InitUsingHamedMethod() */
 
/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PADDS.
******************************************************************************/
void ParaPADDS::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Search for the pareto front using PADDS.
******************************************************************************/
void ParaPADDS::Optimize(void)
{
   ArchiveStruct tmpArch;
   bool bBanner = false;
   StatusStruct pStatus;
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();

   //variables for MPI/parallel 
   int number_of_times_slave_worked;
   bool work_left;
   int nslaves, slaveindex, nxtsid, iLeft, eval;
   const int tag = 0;
   MPI_Status status;
   const int dowork=101;
   const int stopwork=102;
   bool bSynch = SynchReceives();
   double Pn;
   FILE * pPnFile;

   int signal;

   InitFromFile(GetInFileName());

   m_num_dec = pGroup->GetNumParams();
   m_num_objs = m_pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);
   m_seed = GetRandomSeed();
   MPI_Comm_size(MPI_COMM_WORLD, &m_nprocessors);
   MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
      
   //vector for sends and recvs
   m_stest_flat = new double[m_num_dec+m_num_objs];

   #if(PARA_PADDS_DEBUG == 1)
      printf("Preparing to optimize\n");
      printf("m_num_dec = %d\n", m_num_dec);
      printf("m_num_objs = %d\n", m_num_objs);
      printf("m_seed = %d\n", m_seed);
      printf("m_nprocessors = %d\n", m_nprocessors);
      printf("m_rank = %d\n", m_rank);
      fflush(stdout);
   #endif

   WriteSetup(m_pModel, "ParaPADDS - Parallel Pareto Archived Dynamically Dimensioned Search");
   //write banner
   WriteBanner(m_pModel, "gen   ", "trials remaining");
   
   //number of initial solutions
   nslaves = m_nprocessors - 1;
   int its = (int)(0.005 * m_maxiter); 
   if(its < nslaves) its = nslaves;
   if(m_maxiter < nslaves) m_maxiter += nslaves;

   ArchiveStruct * stest = NULL;
   ArchiveStruct * sbest = NULL;

   double * S_min, * S_max; //low and high bounds
   S_min = new double[m_num_dec];
   S_max = new double[m_num_dec];

   for (int i = 0; i < m_num_dec; i++)
   { 
      S_min[i] = pGroup->GetParamPtr(i)->GetLwrBnd();
      S_max[i] = pGroup->GetParamPtr(i)->GetUprBnd();
   }

   #if(PARA_PADDS_DEBUG == 1)
      printf("Preparing to optimize, phase 2\n");
      printf("its = %d\n", its);
      for (int i = 0; i < m_num_dec; i++)
      { 
         printf("i = %d, S_min = %lf,  S_max = %lf\n", i, S_min[i], S_max[i]);
      }
      fflush(stdout);
   #endif

   //master section of code
   int master = 0;
   int num_sent = 0; //messages sent
   int num_rcvd = 0; //messages received
   if (m_rank == master)
   {

      pPnFile = fopen("OstPADDSPn.txt", "w");
      fprintf(pPnFile, "EVAL  Pn\n");
      fclose(pPnFile);
      
      slaveindex = 0;
      nxtsid = 0;
    
      #if(PARA_PADDS_DEBUG == 1)
         printf("Proc %d : Entering master section\n", m_rank);
         printf("Proc %d : nslaves = %d\n", m_rank, nslaves);
         printf("Proc %d : slaveindex = %d\n", m_rank, slaveindex);
         fflush(stdout);
      #endif
  
      //must send 1 work message per iteration plus 1 stop work message per slave
      WriteInnerEval(WRITE_DDS, 0, '.');
      for(num_sent = 0; num_sent <= (m_maxiter + nslaves); num_sent++)
      {
         if(IsQuit() == true){ MPI_Abort(MPI_COMM_WORLD, 0);}

         #if(PARA_PADDS_DEBUG == 1)
         printf("Proc %d : Master loop processing\n", m_rank);
         printf("Proc %d : num_sent = %d\n", m_rank, num_sent);
         printf("Proc %d : num_rcvd = %d\n", m_rank, num_rcvd);
         fflush(stdout);
         #endif

         /* -------------------------------------------------------------------
         If more messages than slaves have been sent, it is time to start 
         receiving the results and processing them.
         ------------------------------------------------------------------- */
         if(num_sent >= nslaves)
         {
            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : Master receive processing\n", m_rank);
               printf("Proc %d : num_sent = %d , nslaves = %d\n", m_rank, num_sent, nslaves);
               fflush(stdout);
            #endif

            if(bSynch == true)
            {
               slaveindex = nxtsid + 1;
               nxtsid = (nxtsid + 1) % nslaves;
            }
            else
            {
               slaveindex = MPI_ANY_SOURCE;
            }            

            MPI_Recv(m_stest_flat, m_num_dec+m_num_objs,MPI_DOUBLE,slaveindex,tag,MPI_COMM_WORLD,&status); 

            //prepare permanent archive entry
            ArchiveStruct * stest = new ArchiveStruct;
            stest->nF = m_num_objs;
            stest->nX = m_num_dec;
            stest->F = new double[m_num_objs];
            stest->X = new double[m_num_dec];
            stest->Z = -999.999;
            stest->pNext = NULL;
            for (int j = 0; j < m_num_dec; j++)
            {
               stest->X[j] = m_stest_flat[j];
            }
            for (int j = 0; j < m_num_objs; j++)
            {
               stest->F[j] = m_stest_flat[j+m_num_dec];
            }      

            //determine source
            slaveindex = status.MPI_SOURCE;                                     

            num_rcvd++;

            WriteInnerEval(num_rcvd, m_maxiter, '.');

            pStatus.curIter = num_rcvd;
            pStatus.maxIter = m_maxiter;
            pStatus.pct = (float)100.00*(float)num_rcvd/(float)m_maxiter;
            pStatus.numRuns = num_rcvd;
        WriteStatus(&pStatus);

            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : Master received the following result from processor %d :\n", m_rank, slaveindex);
               for (int j = 0; j < m_num_dec; j++)
               {
                  printf("Proc %d : X[%d] = %lf\n", m_rank, j, stest->X[j]);
               }
               for (int j = 0; j < m_num_objs; j++)
               {
                  printf("Proc %d : F[%d] = %lf\n", m_rank, j, stest->F[j]);
               }      
               fflush(stdout);
            #endif

            /* -------------------------------------------------------------------
            Determine if we are receiving results from the initialization phase or 
            the search phase and process results accordingly.
            ------------------------------------------------------------------- */
            if(num_rcvd < its) //receiving intitialization results
            {
               #if(PARA_PADDS_DEBUG == 1)
                  printf("Proc %d : Master received initialization results\n", m_rank);
                  printf("Proc %d : num_rcvd = %d , its = %d\n", m_rank, num_rcvd, its);
                  fflush(stdout);
               #endif

               //process result stored in m_stest_flat
               int result = UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);

               #if(PARA_PADDS_DEBUG == 1)
                  printf("Proc %d : Archive updated, result = %d\n", m_rank, result);
                  fflush(stdout);
               #endif

               if(result == ARCHIVE_NON_DOM)
               {
                  #if(PARA_PADDS_DEBUG == 1)
                     printf("Proc %d : Found new non-dominated solution, printing revised list.\n", m_rank);
                     fflush(stdout);
                  #endif
                  WriteMultiObjRecord(m_pModel, num_rcvd, m_pNonDom, (double)(m_maxiter-num_rcvd));
               }
            }/* end if() */
            else //receiving search results
            {
               #if(PARA_PADDS_DEBUG == 1)
                  printf("Proc %d : Master received PADDS search results\n", m_rank);
                  printf("Proc %d : num_rcvd = %d , its = %d\n", m_rank, num_rcvd, its);
                  fflush(stdout);
               #endif

               //checking to see if x_curr dominates x_new
               bool sbest_dominates_stest = false;

               /* --------------------------------------------------
               Could be a mistmatch between the sbest tested here 
               and the one used in preparing the original send.
               -------------------------------------------------- */
               #if(PARA_PADDS_DEBUG == 1)
                  printf("Proc %d : Checking dominion status stest vs. sbest\n", m_rank);
                  for (int j = 0; j < m_num_objs; j++)
                  {
                     printf("Proc %d : j = %d, Ftest = %lf , Fbest = %lf\n", m_rank, j, stest->F[j], sbest->F[j]);
                  }
                  fflush(stdout);
               #endif
               if(dominion_status(stest, sbest) == 2)
               {
                  #if(PARA_PADDS_DEBUG == 1)
                     printf("Proc %d : dominion_status = 2 (sbest_dominates_stest)\n", m_rank);
                     printf("Proc %d : Adding dominated solution to archive\n", m_rank);
                     printf("Proc %d : Marking dominance flag (-1)\n", m_rank);
                     fflush(stdout);
                  #endif
                  sbest_dominates_stest = true;
                  m_dominance_flag = -1;
                  UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);
               }/* new solution not dominated by sbest, but what about others in the archive? */
               else
               {
                  #if(PARA_PADDS_DEBUG == 1)
                     printf("Proc %d : dominion_status != 2 (stest in non-dminated wrt sbest)\n", m_rank);
                     printf("Proc %d : checking for duplicate archive entry\n", m_rank);
                     fflush(stdout);
                  #endif

                  //check to see if it is a duplicate, if yes, flag = 0, do not Update_Archive, and do not add
                  bool stest_is_duplicate = true;
                  for (ArchiveStruct * pCur = m_pNonDom; pCur != NULL; pCur = pCur->pNext)
                  {
                     for(int j = 0; j < pCur->nF; j++)
                     {
                        if (stest->F[j] != pCur->F[j])
                        {
                           stest_is_duplicate = false;
                           break;
                        }
                     }
                     if(stest_is_duplicate == false)
                     {
                        #if(PARA_PADDS_DEBUG == 1)
                           printf("Proc %d : solution IS NOT a duplicate\n", m_rank);
                           fflush(stdout);
                        #endif
                        break;
                     }
                  }/* end for() */

                  if(stest_is_duplicate == true)
                  {
                     #if(PARA_PADDS_DEBUG == 1)
                        printf("Proc %d : solution IS a duplicate\n", m_rank);
                        printf("Proc %d : clearing dominance flag and destroying solution\n", m_rank);
                        fflush(stdout);
                     #endif
                     m_dominance_flag = 0;
                     DestroyArchive(stest);
                  }

                  if (stest_is_duplicate == false)
                  {
                     #if(PARA_PADDS_DEBUG == 1)
                        printf("Proc %d : Adding to archive and checking dominance relation.\n", m_rank);
                        fflush(stdout);
                     #endif

                     int result = UpdateArchive(stest->X, stest->nX, stest->F, stest->nF);
                     if(result == ARCHIVE_NON_DOM)
                     {
                        #if(PARA_PADDS_DEBUG == 1)
                           printf("Proc %d : A new non-domnated solution has been found. Printing revised list.\n", m_rank);
                           fflush(stdout);
                        #endif

                        WriteInnerEval(WRITE_ENDED, 0, '.');
                        WriteMultiObjRecord(m_pModel, (num_rcvd), m_pNonDom, (double)(m_maxiter-num_rcvd));
                        if((m_maxiter-num_rcvd) > 0)
                        {
                           WriteInnerEval(WRITE_DDS, 0, '.');
                        }
                        bBanner = true;
                     }
                     else
                     {
                        #if(PARA_PADDS_DEBUG == 1)
                           printf("Proc %d : The new solution is dominated by existing solutions. Clearing banner print flag.\n", m_rank);
                           fflush(stdout);
                        #endif

                        bBanner = false;
                     }
                  }/* end else() */

                  #if(PARA_PADDS_DEBUG == 1)
                     printf("Proc %d : Clearing duplicate flag and check for Calc_Z(). \n", m_rank);
                     printf("Proc %d : m_dominance_flag = %d\n", m_rank, m_dominance_flag);
                     fflush(stdout);
                  #endif

                  stest_is_duplicate = false;

                  if(m_dominance_flag != -1)
                  {
                     #if(PARA_PADDS_DEBUG == 1)
                        printf("Proc %d : Calculating revised Z metric for non-dominated solutions. \n", m_rank);
                        fflush(stdout);
                     #endif

                     Calc_Z(m_pNonDom);
                  }
               }/* end else(dominion status) */
            }/* end if() */
         }/* end if() */
         /* -------------------------------------------------------------------
         Fewer messages than slaves have been sent, increment slaveindex so that
         subsequent work assignment goes to a new slave.
         ------------------------------------------------------------------- */
         else
         {
            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : Fewer messages (%d) than slaves (%d).\n", m_rank, num_sent, nslaves);
               printf("Proc %d : Incrementing slave index (currently is %d)\n", m_rank, slaveindex);
               fflush(stdout);
            #endif

            slaveindex++;
         }/* end else() */

         /* -------------------------------------------------------------------
         Determine what type of work message to send. Either a new assignment
         or a stop work message.
         ------------------------------------------------------------------- */
         #if(PARA_PADDS_DEBUG == 1)
            printf("Proc %d : Determining the type of work message to send.\n", m_rank);
            printf("Proc %d : num_sent = %d, num_rcvd = %d, m_maxiter = %d\n", m_rank, num_sent, num_rcvd, m_maxiter);
            fflush(stdout);
         #endif
         //have we sent out enough work assignements?
         if(num_sent <= m_maxiter) //more work to be done
         {
            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : More work to be done.\n", m_rank);
               fflush(stdout);
            #endif

            /* -------------------------------------------------------------------
            Determine if we are sending work from the initialization phase or 
            from the search phase and prepare work accordingly.
            ------------------------------------------------------------------- */
            if(num_sent < its) //sending intitialization tasks
            {
               if(m_InitTaskMethod == INIT_TASK_URANDOM)
               {
                  //generate a new random candidate
                  for (int j = 0; j < m_num_dec; j++)
                  {
                     m_stest_flat[j] = S_min[j] + (S_max[j] - S_min[j])*UniformRandom();
                  }/* end for() */
               }
               else if (m_InitTaskMethod == INIT_TASK_HAMED)
               {
                  InitUsingHamedMethod(m_num_dec, m_stest_flat, S_min, S_max, m_HamedOffset, slaveindex, nslaves);
               }
               else //default to URANDOM
               {
                  //generate a new random candidate
                  for (int j = 0; j < m_num_dec; j++)
                  {
                     m_stest_flat[j] = S_min[j] + (S_max[j] - S_min[j])*UniformRandom();
                  }/* end for() */ 
               }
            }/* end if() */
            else //sending search tasks
            {
               //Calculating initial Selection Metric Z:
               if(num_rcvd == its)
                  Calc_Z(m_pNonDom);

               //DDS bases search on how many iterations are left
               iLeft = m_maxiter - its;

               if (m_dominance_flag == -1)
               {
                  sbest = SelectFrom(m_pNonDom);
               }
               else
               {
                  sbest = m_pNonDom; //Archive[Archive.size() - 1]; 
               }

               // %% DDS
               eval = num_rcvd - its;
               if(eval <= 2*nslaves)
               {
                  Pn = 1.00; //have each slave perturb all parameters the first time through
               }
               else
               {
                  Pn = 1.0 - log((double)(eval-2*nslaves)) / log((double)(iLeft-2*nslaves)); // probability each DV selected
               }

               pPnFile = fopen("OstPADDSPn.txt", "a");
               fprintf(pPnFile, "%04d  %E\n", eval, Pn);
               fclose(pPnFile);

               int dvn_count = 0;

               //copy existing parameter values
               for(int j = 0; j < m_num_dec; j++)
               {
                  m_stest_flat[j] = sbest->X[j];
               }
               for(int j = 0; j < m_num_objs; j++)
               {
                  m_stest_flat[j+m_num_dec] = sbest->F[j];
               }

               //perturb a subset of the parameters
               for (int j = 0; j < m_num_dec; j++)
               {
                  if (UniformRandom() < Pn)
                  {
                     dvn_count += 1;
                     double new_value = neigh_value_continuous(sbest->X[j], S_min[j], S_max[j], m_fraction1);
                     m_stest_flat[j] = new_value;// change relevant dec var value in stest
                  }
               }
               if (dvn_count == 0)
               {         
                  int dec_var = (int)(ceil(m_num_dec * UniformRandom() - 1.00));
                  double new_value = neigh_value_continuous(sbest->X[dec_var], S_min[dec_var], S_max[dec_var], m_fraction1);
                  m_stest_flat[dec_var] = new_value;
               }
            }/* end if() */

            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : Sending work to processor %d.\n", m_rank, slaveindex);
               fflush(stdout);
            #endif
   
            //send out new work assignement
            signal = dowork;
            MPI_Send(&signal,1,MPI_INT,slaveindex,tag,MPI_COMM_WORLD); 
            MPI_Send(m_stest_flat,m_num_dec+m_num_objs,MPI_DOUBLE,slaveindex,tag,MPI_COMM_WORLD);
         }/* end if() */
         else //no more work
         {
            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : No more work, sending stop work to processor %d.\n", m_rank, slaveindex);
               fflush(stdout);
            #endif

            signal = stopwork;
            MPI_Send(&signal,1,MPI_INT,slaveindex,tag,MPI_COMM_WORLD);
         }/* else() */
      } /* end for() */
      if(bBanner == false)
      {
         WriteInnerEval(WRITE_ENDED, 0, '.');
      }
   }/* end if(master) */
   /* slave processing --- receive an assignment and process */
   else
   {
      #if(PARA_PADDS_DEBUG == 1)
         printf("Proc %d : Entering slave section\n", m_rank);
         fflush(stdout);
      #endif

      number_of_times_slave_worked = 0;
      work_left = true;

      #if(PARA_PADDS_DEBUG == 1)
         printf("Proc %d : Initialized number_of_times_slave_worked and work_left\n", m_rank);
         fflush(stdout);
      #endif

      while(work_left == true)
      {

         if(IsQuit() == true){ MPI_Abort(MPI_COMM_WORLD, 0);}

         #if(PARA_PADDS_DEBUG == 1)
            printf("Proc %d : work_left = true\n", m_rank);
            printf("Proc %d : receiving work, work_left = true\n", m_rank);
            fflush(stdout);
         #endif

         MPI_Recv(&signal,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status); 
         if(signal == stopwork)
         {
            work_left = false;
            printf("Proc %d : received termination signal\n", m_rank);
         }
         else if(signal == dowork)
         {
            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : receiving a new work assignment\n", m_rank);
               fflush(stdout);
            #endif

            tmpArch.nF = m_num_objs;
            tmpArch.nX = m_num_dec;
            tmpArch.X = &(m_stest_flat[0]);
            tmpArch.F = &(m_stest_flat[m_num_dec]);

            number_of_times_slave_worked = number_of_times_slave_worked + 1;
            MPI_Recv(m_stest_flat,m_num_dec+m_num_objs,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status); 

            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : running model\n", m_rank);
               fflush(stdout);
            #endif

            //run model
            F(&tmpArch);

            #if(PARA_PADDS_DEBUG == 1)
               printf("Proc %d : sending results\n", m_rank);
               fflush(stdout);
            #endif

            //send back results
            MPI_Send(&(m_stest_flat[0]),m_num_dec+m_num_objs,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
         }
         else
         {
            printf("Proc %d : unknown signal, error\n", m_rank);
            break;
         }
      }/* end while() */

      printf("Proc %d : handled %d tasks\n", m_rank, number_of_times_slave_worked);
   }/* end else() */

   #if(PARA_PADDS_DEBUG == 1)
      printf("Proc %d : All done with main processing loop\n", m_rank);
      fflush(stdout);
   #endif

   if(m_rank == master)
   {
      WriteMultiObjOptimal(m_pModel, m_pNonDom, m_pDom);

      pStatus.pct = 100.00;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //write algorithm metrics
      WriteAlgMetrics(this);

      //reporting the results
      FILE * myfile = fopen("c_results.txt", "w");
      for (ArchiveStruct * pCur = m_pNonDom; pCur != NULL; pCur = pCur->pNext)
      {
         for (int j = 0; j < m_num_objs - 1; j++)
         {
            fprintf(myfile, "%E\t", pCur->F[j]);
         }
         fprintf(myfile, "%E\n", pCur->F[m_num_objs - 1]);
      }
      fclose(myfile); 
   }/* end if(final reporting) */
}/* end Optimize() */

/******************************************************************************
Calc_Z()
******************************************************************************/
void ParaPADDS::Calc_Z(ArchiveStruct * archive)
{
   ArchiveStruct * pA;
   ArchiveStruct ** pSorted;
   int i, archive_size;

   //compute size of archive
   archive_size = 0;
   for(pA = archive; pA != NULL; pA = pA->pNext)
   {   
      archive_size++;
   }
   //make space for sorted archive
   pSorted = new ArchiveStruct *[archive_size];
   //initialize array
   i = 0;
   for(pA = archive; pA != NULL; pA = pA->pNext)
   {
      pSorted[i] = pA;
      i++;   
   }   

   switch (m_Select_metric)
   {
      case 0://RND
      {
         for(pA = archive; pA != NULL; pA = pA->pNext)
         {   
            pA->Z = 1.0;
         }
         break;
      }
      case 1://CD
      {
         for(pA = archive; pA != NULL; pA = pA->pNext)
         {
            pA->Z = 0.0;
         }
         for (int obj = 0; obj < archive->nF; obj++)
         {
            m_dim = obj;

            SortArchive(pSorted, archive_size, obj);

            for (int i = 1; i < archive_size - 1; i++)
            {
               double F0 = pSorted[0]->F[obj];
               double F1 = pSorted[i - 1]->F[obj];
               double F2 = pSorted[i + 1]->F[obj];
               double F4 = pSorted[archive_size - 1]->F[obj];
               pSorted[i]->Z += fabs(F1 - F2) / fabs(F0 - F4);
            }

            if(archive_size > 1)
              pSorted[0]->Z = pSorted[1]->Z;
            if(archive_size > 2)
              pSorted[archive_size - 1]->Z = pSorted[archive_size - 2]->Z;
         }/* end for(each objective function) */
         break;
      }/* end case() */
      case 2://HVC_ESTIMATE
      {
         for (int i = 0; i < archive_size; i++)
         {
            pSorted[i]->Z = 0;
         }

         //define boundries
         double * f_low_bound = new double[m_num_objs];
         double * f_high_bound = new double[m_num_objs];

         //compute boundaries
         for (int i = 0; i < m_num_objs; i++)
         {
            m_dim = i;
            SortArchive(pSorted, archive_size, i);

            f_low_bound[i] = pSorted[0]->F[i];
            f_high_bound[i] = pSorted[archive_size - 1]->F[i];
         }

         //prepare monte carlo archive
         int dots_num = 100;
         ArchiveStruct ** mc_points = new ArchiveStruct*[dots_num];

         for (int i = 0; i < dots_num; i++)
         {
            //allocate sample
            ArchiveStruct * dot = new ArchiveStruct;
            dot->F = new double[m_num_objs];
            dot->nF = m_num_objs;
            dot->X = new double[m_num_dec];
            dot->nX = m_num_dec;
            dot->pNext = NULL;
            dot->Z = -999.999;
            
            for (int j = 0; j < m_num_objs; j++)
            {                                 
               dot->F[j] = f_low_bound[j] + (f_high_bound[j] - f_low_bound[j])*UniformRandom();
            }

            mc_points[i] = dot;
         }/* end for() */

         for (int i = 0; i < dots_num; i++)
         {
            int jj = archive_size; 
            bool any_good = false;

            for (int j = 0; j < archive_size; j++)
            {
               if (dominion_status(pSorted[j], mc_points[i]) == 1)
               {
                  jj = j;
                  any_good = true;
                  break;
               }
            }

            if (any_good == true)
            {
               for (int k = jj + 1; k < archive_size; k++)
               {
                  if (dominion_status(pSorted[k], mc_points[i]) == 1)
                  {
                     //this dot means nothing, goto the next dot.
                     any_good = false;
                     break;
                  }
               }
            }

            //this dot is only dominated by archive[jj]
            if (any_good == true)
            {
               pSorted[jj]->Z += 1.00;
               any_good = false;
            }
         }/* end for() */

         //normalize z
         double best_z = 0.00;

         for (int i = 0; i < archive_size; i++)
         {
            pSorted[i]->Z = (pSorted[i]->Z / (double)dots_num);
            if (pSorted[i]->Z > best_z)
            {
               best_z = pSorted[i]->Z;
            }
         }/* end for() */

         for (int i = 0; i < archive_size; i++)
         {
            if (pSorted[i]->Z == 0)
            {
               pSorted[i]->Z = 0.5*best_z;
            }
         }

         delete [] f_low_bound;
         delete [] f_high_bound;
         for(int i = 0; i < dots_num; i++)
         {
            DestroyArchive(mc_points[i]);
         }
         delete [] mc_points;
         break;
      }/* end case() */
      case 3://HVC_EXACT
      {
         int dataNumber = archive_size;
         int dimension = m_num_objs;

         double * refPoint = new double[dimension];

         for(int i = 0; i < dimension; i++)
         {
            m_dim = i;
            SortArchive(pSorted, archive_size, i);
            refPoint[i] = 1.00001*(pSorted[archive_size - 1]->F[i]);
         }/* end for() */

         double ** pointsInitial = new double *[dataNumber];
         for (int i = 0; i < dataNumber; i++)
         {
            pointsInitial[i] = new double[dimension];

            for (int j = 0; j < m_num_objs; j++)
            {
               pointsInitial[i][j] = pSorted[i]->F[j];
            }
         }

         double HyperVolume = HV(dataNumber, dimension, refPoint, pointsInitial);

         //destroy pointsInitial
         for (int j = 0; j < dataNumber; j++)
         {
            delete[] pointsInitial[j];
         }
         delete[] pointsInitial;

         double best_z = 0;

         for(int i = 0; i < dataNumber; i++) 
         {
            int * included_points = new int[dataNumber];
            for(int j = 0; j < dataNumber; j++)
            {
               included_points[j] = j;
               //leave out ith point
               if(j >= i) included_points[j]++;
            }

            double ** pointsInitial_sub = new double * [dataNumber - 1];
            for (int j = 0; j < dataNumber - 1; j++)
            {
               pointsInitial_sub[j] = new double[dimension];

               for (int k = 0; k < m_num_objs; k++)
               {
                  pointsInitial_sub[j][k] = pSorted[included_points[j]]->F[k];
               }
            }/* end for() */

            pSorted[i]->Z = HyperVolume - HV(dataNumber - 1, dimension, refPoint, pointsInitial_sub);

            if(pSorted[i]->Z > best_z)
            {
               best_z = pSorted[i]->Z; 
            }

            //destroy pointsInitial_sub
            for (int j = 0; j < dataNumber - 1; j++)
            {
               delete[] pointsInitial_sub[j];
            }
            delete[] pointsInitial_sub;
            delete [] included_points;
         }/* end for() */

         // taking care of the edges
         for (int i = 0; i < dimension; i++)
         {
            m_dim = i;
            SortArchive(pSorted, archive_size, i);

            pSorted[0]->Z = best_z;
            pSorted[archive_size - 1]->Z = best_z;
         }

         delete [] refPoint;
         break;
      }/* end case() */
   }/* end switch() */

   delete [] pSorted;
}/* end Calc_Z() */

/******************************************************************************
SortPoints()

Sort the array based on the which index.
******************************************************************************/
void ParaPADDS::SortPoints(double ** X, int size, int which)
{
   double F1, F2;
   double * pTmp;
   for(int i = 0; i < size; i++)
   {
      for(int j = (i+1); j < size; j++)
      {
         F1 = X[i][which];
         F2 = X[j][which];
         if(F2 < F1)
         {
            pTmp = X[i];
            X[i] = X[j];
            X[j] = pTmp;
         }
      }
   }/* end for() */
}/* end SortPoints() */

/******************************************************************************
bool_vec_to_ulong()

Convert bitwise array of bools into equivalent integer.

to_ulong() simply calculates bit[0]*2^0 + bit[1]*2^1 + bit[2]*2^2
******************************************************************************/
int ParaPADDS::bool_vec_to_ulong(bool * pB, int size)
{
   int sum = 0;
   for(int i = 0; i < size; i++)
   {
      sum += (int)(pB[i]*pow((double)2, i));
   }/* end for() */
   return sum;
}/* end bool_vec_to_ulong() */

/******************************************************************************
ulong_to_bool_vec()

Convert integer into equivalent bitwise array of bools.
******************************************************************************/
void ParaPADDS::ulong_to_bool_vec(int val, bool * pB, int size)
{
   int mask = 1;
   int tmp;
   for(int i = 0; i < size; i++)
   {
      tmp = (mask & val);
      if(tmp == 0)
        pB[i] = false;
      else
        pB[i] = true;
      mask = mask << 1;
   }/* end for() */
}/* end ulong_to_bool_vec() */

/******************************************************************************
SortArchive()

Sort the archive by objective function or by Z if whichObj == -1.
******************************************************************************/
void ParaPADDS::SortArchive(ArchiveStruct ** pArch, int size, int whichObj)
{
   double F1, F2;
   ArchiveStruct * pTmp;
   for(int i = 0; i < size; i++)
   {
      for(int j = (i+1); j < size; j++)
      {
         if(whichObj == -1)
         {
            F1 = pArch[i]->Z;
            F2 = pArch[j]->Z;
         }
         else
         {
            F1 = pArch[i]->F[whichObj];
            F2 = pArch[j]->F[whichObj];
         }
         if(F2 < F1)
         {
            pTmp = pArch[i];
            pArch[i] = pArch[j];
            pArch[j] = pTmp;
         }
      }
   }/* end for() */
}/* end SortArchive() */

/******************************************************************************
DestroyArchive()

Free up memory of archive.
******************************************************************************/
void ParaPADDS::DestroyArchive(ArchiveStruct * pArch)
{
   ArchiveStruct * pDel;
   for(ArchiveStruct * pCur = pArch; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }
}/* end DestroyArchive() */

/******************************************************************************
UpdateArchive()

Update the dominated and non-dominated archives with latest sample.
******************************************************************************/
int ParaPADDS::UpdateArchive(double * pX, int nX, double * pF, int nF)
{
   int i;
   double Fcur, Ftst;
   ArchiveStruct * pArch, * pCur, * pPrev, * pNxt;
   bool bDominates, bIsDominated, bMarkForInsertion;
   pArch =  new ArchiveStruct;
   pArch->F = pF;
   pArch->X = pX;
   pArch->nX = nX;
   pArch->nF = nF;
   pArch->pNext = NULL;
   pArch->Z = -999.999;
   
   //first entry is always non-dominated
   if((m_NumDom == 0) && (m_NumNonDom == 0))
   {
      m_pDom = NULL;
      m_pNonDom = pArch;
      m_NumNonDom++;
      return ARCHIVE_NON_DOM;
   }

   //assume solution is non-dominated until we discover otherwise
   bMarkForInsertion = true;

   //compare against current list of non-dominated solutions
   ArchiveStruct * pDummy;
   for(pCur = m_pNonDom; pCur != NULL;)
   {
      //save next item since pCur->pNext may be changed during processing
      pDummy = pCur->pNext;

      //does new solution (Ftst) dominate the existing solution (Fcur)?
      bDominates = true;
      for(i = 0; i < pArch->nF; i++)
      {
         Fcur = pCur->F[i];
         Ftst = pArch->F[i];
         if(Fcur < Ftst)
         {
            bDominates = false;
            break;
         }/* end if() */
      }/* end for() */

      //is new solution (Ftst) dominated by an existing solution (Fcur)?
      if(bDominates == false)
      {
         bIsDominated = true;
         for(i = 0; i < pArch->nF; i++)
         {
            Fcur = pCur->F[i];
            Ftst = pArch->F[i];
            if(Ftst < Fcur)
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
      if(bDominates == true)
      {
         //solution to be removed is at head of list.
         if(pCur == m_pNonDom)
         {
            m_pNonDom = pCur->pNext;
            m_NumNonDom--;
         }/* end if() */
         else //somewhere in middle or at end.
         {
            for(pPrev = m_pNonDom; pPrev->pNext != pCur;)
            {
               pPrev = pPrev->pNext;
            }
            pNxt = pCur->pNext;
            pPrev->pNext = pNxt; //removes item from list
            m_NumNonDom--;
         }/* end else() */

         //insert at head of dominated list
         pCur->pNext = NULL;
         pNxt = m_pDom;
         m_pDom = pCur;
         m_pDom->pNext = pNxt;
         m_NumDom++;         
      }/* end if() */
      /* -----------------------------------------------------------------------------
      New solution is dominated. Make note so that it is not inserted into the list
      of non-dominated solutions.
      ----------------------------------------------------------------------------- */
      else if(bIsDominated == true)
      {
         bMarkForInsertion = false;   
      }/* end if() */

      //advance to next item
      pCur = pDummy;
   }/* end for() */

   //insert new solution into list of non-dominated solutions?
   if(bMarkForInsertion == true)
   {
      //insert at head of non-dominated list
      pNxt = m_pNonDom;
      m_pNonDom = pArch;
      m_pNonDom->pNext = pNxt;
      m_NumNonDom++;               
      return ARCHIVE_NON_DOM;
   }/* end if() */
   else
   {
      //insert at head of dominated list
      pNxt = m_pDom;
      m_pDom = pArch;
      m_pDom->pNext = pNxt;
      m_NumDom++;               
      return ARCHIVE_DOM;
   }/* end else() */
}/* end UpdateArchive() */

/******************************************************************************
dominion_status()

Determine whether solution x1 dominates solution x2.
******************************************************************************/
int ParaPADDS::dominion_status(ArchiveStruct * x1, ArchiveStruct * x2)
{
   int ds = 1;

   for (int i = 0; i < m_num_objs; i++)
   {
      if (x1->F[i] > x2->F[i])
      {
         goto try_2;
      }
   }
   return 1;

try_2:

   ds = 2;
   for (int i = 0; i < m_num_objs; i++)
   {
      if (x1->F[i] < x2->F[i])
      {
         return 0;
      }
   }

   return 2;
}/* end dominion_status() */

/******************************************************************************
SelectFrom()

Select and entry from the archive.
******************************************************************************/
ArchiveStruct * ParaPADDS::SelectFrom(ArchiveStruct * pArchive)
{
   ArchiveStruct * pA;
   ArchiveStruct ** archive;
   int i, archive_size;

   //compute size of archive
   archive_size = 0;
   for(pA = pArchive; pA != NULL; pA = pA->pNext)
   {   
      archive_size++;
   }
   //make space for sorted archive
   archive = new ArchiveStruct *[archive_size];
   //initialize array
   i = 0;
   for(pA = pArchive; pA != NULL; pA = pA->pNext)
   {
      archive[i] = pA;
      i++;   
   }   

   double * z_cum = new double[archive_size];
   z_cum[0] = 0;
   for (int i = 0; i < archive_size; i++)
   {
      if (i == 0)
      {
         z_cum[i] = archive[i]->Z;
      }
      else
      {
         z_cum[i] = z_cum[i - 1] + archive[i]->Z;
      }
   }/* end for() */

   double t = UniformRandom() * z_cum[archive_size - 1];

   int ii = 0;

   for (int i = 0; i < archive_size; i++)
   {
      if (z_cum[i] >= t)
      {
         ii = i;
         break;
      }
   }

   delete [] z_cum;
   pA = archive[ii];
   delete [] archive;
   return pA;
}/* end SelectFrom() */

/******************************************************************************
neigh_value_continuous()
******************************************************************************/
double ParaPADDS::neigh_value_continuous(double s, double s_min, double s_max, double r)
{
   double s_range = s_max - s_min;

   double snew = s + GaussRandom() * r * s_range;

   double P_Abs_or_Ref = UniformRandom();

   if (snew < s_min) 
   {
      if (P_Abs_or_Ref <= 0.5)
      {
         snew = s_min + (s_min - snew);
      }
      else
      {
         snew = s_min;
      }
      if (snew > s_max)
      { 
         snew = s_min; 
      }
   }
   else if (snew > s_max)
   {
      if (P_Abs_or_Ref <= 0.5)
      {
         snew = s_max - (snew - s_max);
      }
      else
      {
         snew = s_max;
      }

      if (snew < s_min)
      {
         snew = s_max;
      }
   }
   return snew;
}/* end neigh_value_continuous() */

/******************************************************************************
HV()

Hypervolume calculation.
******************************************************************************/
double ParaPADDS::HV(int data_n, int dim_n, double * ref, double ** points)
{
   int i, j;
   int dataNumber = data_n;

   //set global
   m_dimension = dim_n;

   double ** pointsInitial = new double *[dataNumber];
   for (int n = 0; n < dataNumber; n++)
   {
      pointsInitial[n] = new double[dim_n];
      for (int i = 0; i < dim_n; i++)
      {
         pointsInitial[n][i] = points[n][i];
      }
   }
   double* refPoint = new double[dim_n];
   for (i = 0; i < dim_n; i++)
   {
      refPoint[i] = ref[i];
   }

   // initialize volume
   m_volume = 0.0;
   // sqrt of dataNumber
   m_dSqrtDataNumber = sqrt((double)dataNumber);

   // initialize region
   double* regionLow = new double[dim_n - 1];
   double* regionUp = new double[dim_n - 1];
   for (j = 0; j < dim_n - 1; j++) 
   {
      // determine minimal j coordinate
      double min = NEARLY_HUGE;
      for (i = 0; i < dataNumber; i++)
      {
         if (pointsInitial[i][j] < min) 
         {
            min = pointsInitial[i][j];
         }
      }
      regionLow[j] = min;
      regionUp[j] = refPoint[j];
   }

   // sort pointList according to d-th dimension
   SortPoints(pointsInitial, dataNumber, dim_n - 1);

   // call stream initially
   stream(regionLow, regionUp, pointsInitial, dataNumber, 0, refPoint[m_dimension - 1]);

   for (int n = 0; n < dataNumber; n++)
   {
      delete [] pointsInitial[n];
   }
   delete [] pointsInitial;
   delete [] refPoint;
   delete [] regionLow;
   delete [] regionUp;

   // print hypervolume
   return m_volume;
} /* end HV() */

/******************************************************************************
F()

The objective functions.
******************************************************************************/
void ParaPADDS::F(ArchiveStruct * pA)
{
   m_pModel->GetParamGroupPtr()->WriteParams(pA->X);
   m_pModel->Execute(pA->F, pA->nF);
}/* end F() */

/******************************************************************************
covers()
******************************************************************************/
bool ParaPADDS::covers(double* cub, double * regLow) 
{
   int i_hvc;
   for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
      if (cub[i_hvc] > regLow[i_hvc]) 
      {
         return false;
      }
   }
   return true;
}/* end covers() */

/******************************************************************************
partCovers()
******************************************************************************/
bool ParaPADDS::partCovers(double* cub, double * regUp)
{
   int i_hvc;
   for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++)
   {
      if (cub[i_hvc] >= regUp[i_hvc])
      {
         return false;
      }
   }
   return true;
}/* end partCovers() */

/******************************************************************************
containsBoundary()
******************************************************************************/
int ParaPADDS::containsBoundary(double* cub, double * regLow, int split) 
{
   // condition only checked for split>0
   if (regLow[split] >= cub[split])
   {
      // boundary in dimension split not contained in region, thus
      // boundary is no candidate for the splitting line
      return -1;
   }
   else 
   {
      int j_hvc;
      for (j_hvc = 0; j_hvc < split; j_hvc++) 
      { // check boundaries
         if (regLow[j_hvc] < cub[j_hvc]) 
         {
            // boundary contained in region
            return 1;
         }
      }
   }
   // no boundary contained in region
   return 0;
}/* end containsBoundary() */

/******************************************************************************
getMeasure()
******************************************************************************/
double ParaPADDS::getMeasure(double * regLow, double * regUp) 
{
   double vol_hvc;
   int i_hvc;
   vol_hvc = 1.0;
   for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
      vol_hvc *= (regUp[i_hvc] - regLow[i_hvc]);
   }
   return vol_hvc;
}/* end getMeasure() */

/******************************************************************************
isPile()
******************************************************************************/
int ParaPADDS::isPile(double* cub, double * regLow, double * regUp) 
{
   int pile_hvc;
   int k_hvc;

   pile_hvc = m_dimension;
   // check all dimensions of the node
   for (k_hvc = 0; k_hvc < m_dimension - 1; k_hvc++) 
   {
      // k-boundary of the node's region contained in the cuboid? 
      if (cub[k_hvc] > regLow[k_hvc]) 
      {
         if (pile_hvc != m_dimension) 
         {
            // second dimension occured that is not completely covered
            // ==> cuboid is no pile
            return -1;
         }
         pile_hvc = k_hvc;
      }
   }
   // if pile == this.dimension then
   // cuboid completely covers region
   // case is not possible since covering cuboids have been removed before

   // region in only one dimenison not completly covered 
   // ==> cuboid is a pile 
   return pile_hvc;
}/* isPile() */

/******************************************************************************
computeTrellis()
******************************************************************************/
double ParaPADDS::computeTrellis(double * regLow, double * regUp, double * trellis) 
{
   int i_hvc, j_hvc;
   double vol_hvc;
   int numberSummands_hvc;
   double summand_hvc;
   bool * bitvector_hvc = new bool[m_dimension];

   vol_hvc = 0.0;
   summand_hvc = 0.0;
   numberSummands_hvc = 0;

   // calculate number of summands
   bool * nSummands = new bool[m_dimension];
   for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
      nSummands[i_hvc] = 1;
   }
   numberSummands_hvc = bool_vec_to_ulong(nSummands, m_dimension);

   double* valueTrellis = new double[m_dimension - 1];
   double* valueRegion = new double[m_dimension - 1];
   for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
      valueTrellis[i_hvc] = trellis[i_hvc] - regUp[i_hvc];
   }
   for (i_hvc = 0; i_hvc < m_dimension - 1; i_hvc++) 
   {
      valueRegion[i_hvc] = regUp[i_hvc] - regLow[i_hvc];
   }

   double* dTemp = new double[numberSummands_hvc / 2 + 1];

   // sum
   for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) 
   {

      // set bitvector length to fixed value 16
      // TODO Warning: dimension-1 <= 16 is assumed
      ulong_to_bool_vec(i_hvc, bitvector_hvc, m_dimension);

      // construct summand
      // 0: take factor from region
      // 1: take factor from cuboid
      summand_hvc = 1.0;
      for (j_hvc = 0; j_hvc < m_dimension - 2; j_hvc++) 
      {
         if (bitvector_hvc[j_hvc]) 
         {
            summand_hvc *= valueTrellis[j_hvc];
         }
         else 
         {
            summand_hvc *= valueRegion[j_hvc];
         }
      }
      summand_hvc *= valueRegion[m_dimension - 2];

      // determine sign of summand
      vol_hvc -= summand_hvc;
      dTemp[i_hvc] = -summand_hvc;

      // add summand to sum
      // sign = (int) pow((double)-1, (double)counterOnes+1);
      // vol += (sign * summand); 
   }

   ulong_to_bool_vec(i_hvc, bitvector_hvc, m_dimension);

   summand_hvc = 1.0;
   for (j_hvc = 0; j_hvc < m_dimension - 1; j_hvc++) 
   {
      if (bitvector_hvc[j_hvc]) 
      {
         summand_hvc *= valueTrellis[j_hvc];
      }
      else 
      {
         summand_hvc *= valueRegion[j_hvc];
      }
   }
   vol_hvc -= summand_hvc;

   for (i_hvc = 1; i_hvc <= numberSummands_hvc / 2; i_hvc++) 
   {
      summand_hvc = dTemp[i_hvc];
      summand_hvc *= regUp[m_dimension - 2] - trellis[m_dimension - 2];
      summand_hvc /= valueRegion[m_dimension - 2];
      vol_hvc -= summand_hvc;
   }
   //she has already set it{
   delete [] valueTrellis;
   delete [] valueRegion;
   delete [] dTemp;
   delete [] bitvector_hvc;
   delete [] nSummands;
   //she has already set it}

   return vol_hvc;
}/* end trellis() */

/******************************************************************************
getMedian()

return median of the list of boundaries considered as a set
******************************************************************************/
double ParaPADDS::getMedian(double * bounds, int size) 
{
   double * pSort;
   double tmp, median;

   // do not filter duplicates   
   if (size == 1) 
   {
      return bounds[0];
   }
   else if (size == 2) 
   {
      return bounds[1];
   }

   pSort = new double[size];
   for(int i = 0; i < size; i++)
   {
      pSort[i] = bounds[i];
   }
   for(int i = 0; i < size; i++)
   {
      for(int j = i+1; j < size; j++)
      {
         if(pSort[j] < pSort[i])
         {
            tmp = pSort[j];
            pSort[j] = pSort[i];
            pSort[i] = tmp;
         }
      }
   }

   median = pSort[size/2];

   delete [] pSort;

   return median;
}/* end getMedian() */

/******************************************************************************
stream()

recursive calculation of hypervolume
******************************************************************************/
void ParaPADDS::stream(double * regionLow, double * regionUp, double ** points, int npoints, 
            int split, double cover) 
{
   //--- init --------------------------------------------------------------//
   double coverOld_hvc;
   coverOld_hvc = cover;
   unsigned int coverIndex_hvc = 0;
   int c_hvc;

   //--- cover -------------------------------------------------------------//

   // identify first covering cuboid
   double dMeasure = getMeasure(regionLow, regionUp);
   while (cover == coverOld_hvc && coverIndex_hvc < (unsigned int)npoints)
   {
      if (covers(points[coverIndex_hvc], regionLow)) 
      {
         // new cover value
         cover = points[coverIndex_hvc][m_dimension - 1];
         m_volume += dMeasure * (coverOld_hvc - cover);
      }
      else coverIndex_hvc++;
   }/* end while() */

   /* coverIndex shall be the index of the first point in points which
   * is ignored in the remaining process
   *
   * It may occur that that some points in front of coverIndex have the same
   * d-th coordinate as the point at coverIndex. This points must be discarded
   * and therefore the following for-loop checks for this points and reduces
   * coverIndex if necessary.
   */
   for (c_hvc = coverIndex_hvc; c_hvc > 0; c_hvc--) 
   {
      if (points[c_hvc - 1][m_dimension - 1] == cover) 
      {
         coverIndex_hvc--;
      }
   }/* end for() */

   // abort if points is empty
   if (coverIndex_hvc == 0) 
   {
      return;
   }
   // Note: in the remainder points is only considered to index coverIndex

   //--- allPiles  ---------------------------------------------------------//

   bool allPiles = true;
   unsigned int iii;

   int* piles_hvc = new int[coverIndex_hvc];
   for (iii = 0; iii < coverIndex_hvc; iii++) 
   {
      piles_hvc[iii] = isPile(points[iii], regionLow, regionUp);
      if (piles_hvc[iii] == -1) 
      {
         allPiles = false;

         //she has already set it{
         delete[] piles_hvc;
         //she has already set it{

         break;
      }
   }/* end for() */

   /*
   * tre llis[i] contains the values of the minimal i-coordinate of
   * the i-piles.
   * If there is no i-pile the default value is the upper bpund of the region.
   * The 1-dimensional KMP of the i-piles is: reg[1][i] - tre llis[i]
   *
   */
   if (allPiles) 
   { // sweep

      // initialize trellis with region's upper bound
      double* trellis = new double[m_dimension - 1];
      for (c_hvc = 0; c_hvc < m_dimension - 1; c_hvc++) 
      {
         trellis[c_hvc] = regionUp[c_hvc];
      }

      double current = 0.0;
      double next = 0.0;
      iii = 0;
      do 
      { // while(next != coverNew)
         current = points[iii][m_dimension - 1];
         do 
         { // while(next == current)
            if (points[iii][piles_hvc[iii]] < trellis[piles_hvc[iii]]) 
            {
               trellis[piles_hvc[iii]] = points[iii][piles_hvc[iii]];
            }
            iii++; // index of next point
            if (iii < coverIndex_hvc) 
            {
               next = points[iii][m_dimension - 1];
            }
            else 
            {
               next = cover;
            }

         } while (next == current);
         m_volume += computeTrellis(regionLow, regionUp, trellis) * (next - current);
      } while (next != cover);
   }
   //--- split -------------------------------------------------------------//
   // inner node of partition tree
   else
   {
      double bound = -1.0;
      int boundaries_size = 0;
      int noBoundaries_size = 0;
      double * boundaries = new double[coverIndex_hvc];
      double * noBoundaries = new double[coverIndex_hvc];

      do 
      {
         for (iii = 0; iii < coverIndex_hvc; iii++) 
         {
            int contained = containsBoundary(points[iii], regionLow, split);
            if (contained == 1) 
            {
               boundaries[boundaries_size] = points[iii][split];
               boundaries_size++;
            }
            else if (contained == 0) 
            {
               noBoundaries[noBoundaries_size] = points[iii][split];
               noBoundaries_size++;
            }
         }

         if (boundaries_size >  0) 
         {
            bound = getMedian(boundaries, boundaries_size);
         }
         else if (noBoundaries_size > m_dSqrtDataNumber) 
         {
            bound = getMedian(noBoundaries, noBoundaries_size);
         }
         else 
         {
            split++;
         }
      } while (bound == -1.0);

      //destroy boundaries vector
      delete [] boundaries;
      delete [] noBoundaries;

      double dLast;
      int pointsChild_size = 0;
      double ** pointsChild = new double *[coverIndex_hvc];

      // left child
      // reduce maxPoint
      dLast = regionUp[split];
      regionUp[split] = bound;
      for (iii = 0; iii < coverIndex_hvc; iii++) 
      {
         if (partCovers(points[iii], regionUp)) 
         {
            pointsChild[pointsChild_size] = points[iii];
            pointsChild_size++;
         }
      }
      if (pointsChild_size > 0) 
      {
         stream(regionLow, regionUp, pointsChild, pointsChild_size, split, cover);
      }

      // right child
      // increase minPoint
      pointsChild_size = 0;

      regionUp[split] = dLast;
      dLast = regionLow[split];
      regionLow[split] = bound;
      for (iii = 0; iii < coverIndex_hvc; iii++) 
      {
         if(partCovers(points[iii], regionUp)) 
         {
            pointsChild[pointsChild_size] = points[iii];
            pointsChild_size++;
         }
      }
      if (pointsChild_size > 0) 
      {
         stream(regionLow, regionUp, pointsChild, pointsChild_size, split, cover);
      }
      regionLow[split] = dLast;

      delete [] pointsChild;
   }// end inner node
} /* end stream() */

/******************************************************************************
PARA_PADDS_Program()

Calibrate or optimize the model using PADDS.
******************************************************************************/
void PARA_PADDS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("ParaPADDS", 1);
   ParaPADDS * TheAlg = new ParaPADDS(model);
   MEM_CHECK(TheAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end PARA_PADDS_Program() */
