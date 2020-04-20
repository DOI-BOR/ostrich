/******************************************************************************
File     : SurrogateDbase.cpp
Author   : L. Shawn Matott
Copyright: 2006, L. Shawn Matott

Manages a database of model runs.

Version History
04-18-06    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <stdio.h>
#include <math.h>

#include "SurrogateDbase.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "Utility.h"

/*****************************************************************************
CTOR
*****************************************************************************/
SurrogateDbase::SurrogateDbase(int size, int psize, int n_models)
{
   int i, j;

   m_CurSize = 0;
   m_MaxSize = size;
   m_NumParams = psize;
   m_NumModels = n_models;
   m_pDbase = new DbaseEntry[size];
   m_pAvgRunTimes = new double[n_models];
   MEM_CHECK(m_pAvgRunTimes);

   for(i = 0; i < size; i++)
   {
      m_pDbase[i].run_time = 0.00;
      m_pDbase[i].id = -1;
      m_pDbase[i].F = NEARLY_HUGE;
      m_pDbase[i].pParams = new double[psize];      

      for(j = 0; j  < psize; j++)
      {
         m_pDbase[i].pParams[j] = 0.00;
      }
   }/* end for() */

   m_Temp.pParams = new double[psize];
   for(j = 0; j  < psize; j++)
   {
      m_Temp.pParams[j] = 0.00;
   }

   IncCtorCount();
}/* end CTOR */

/*****************************************************************************
Destroy()
*****************************************************************************/
void SurrogateDbase::Destroy(void)
{
   int i;

   for(i = 0; i < m_MaxSize; i++)
   {
      delete [] (m_pDbase[i].pParams);
   }
   delete [] m_pDbase;

   delete [] m_Temp.pParams;
   delete [] m_pAvgRunTimes;

   IncDtorCount();
}/* end Destroy() */

/* ****************************************************************************
GetNearestNeighbor()

Retrieve the WSSE value of the nearest neighboring entry in the database for 
the given model id.
*****************************************************************************/
double SurrogateDbase::GetNearestNeighbor(int id, double * pX)
{
  int i,j, max_idx;
  double D, Dnear = -1.00;
  double F, Fnear;
  double x1, x2;
  bool first = true;

  max_idx = m_CurSize;
  if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

  for(j = 0; j < max_idx; j++)
    {
      if(m_pDbase[j].id == id)
	{
	  F = m_pDbase[j].F;

	  D = 0.00;
	  for(i = 0; i < m_NumParams; i++){
            x1 = pX[i];
            x2 = m_pDbase[j].pParams[i];
            D += ((x2-x1)*(x2-x1));}
	  D = sqrt(D);

	  if(first == true)
	    {
	      Fnear = F;
	      Dnear = D;
	      first = false;
	    }
	  else
	    {
	      if(D < Dnear)
		{
		  Fnear = F;
		  Dnear = D;
		}
	    }
	}/* end if() */
    }/* end for() */

  return Fnear;
}/* end GetNearestNeighbor() */

/* ****************************************************************************
InvDistWSSE()

Compute an interpolated WSSE value using inverse distance weighting.
*****************************************************************************/
double SurrogateDbase::InvDistWSSE(int id, double * pX)
{
   int i,j, max_idx;
   double Di, Dtot, Dflt, Dmin, Fi, Fest, Wi;
   double x1, x2;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   /* ------------------------------------------
   First compute the total inverse distance
   (Dtot) from the new point (pX) to all other
   points.
   ------------------------------------------ */   
   Dtot = 0.00;
   for(j = 0; j < max_idx; j++)
   {
      if(m_pDbase[j].id == id)
	   {
         Fi = m_pDbase[j].F;
         Di = 0.00;
         for(i = 0; i < m_NumParams; i++){
            x1 = pX[i];
            x2 = m_pDbase[j].pParams[i];
            Di += ((x2-x1)*(x2-x1));}
	      Di = sqrt(Di);

         //point already stored?
         if(Di <= NEARLY_ZERO) return Fi;
         
         Dtot += (1/Di);
      }/* end if() */
   }/* end for() */

   /* ------------------------------------------
   Now compute a filterted Dtot by filtering 
   points with weights that are less than some 
   minimum value (i.e. ignore points that are 
   relatively far away).
   ------------------------------------------ */   
   Dflt = 0.00;
   Dmin = 0.10*Dtot;
   for(j = 0; j < max_idx; j++)
   {
      if(m_pDbase[j].id == id)
	   {   
         Di = 0.00;
         for(i = 0; i < m_NumParams; i++){
            x1 = pX[i];
            x2 = m_pDbase[j].pParams[i];
            Di += ((x2-x1)*(x2-x1));}
	      Di = sqrt(Di);
         Di = 1.00/Di;
         if(Di < Dmin) Di = 0.00;
         Dflt += Di;
      }/* end if() */
   }/* end for() */

   /* -----------------------------------------
   If all points are relatively far away, don't
   filter any of them.
   ----------------------------------------- */
   if(Dflt <= NEARLY_ZERO){
      Dflt = Dtot;
      Dmin = 0.00;}

   /* ------------------------------------------
   Compute filtered weights and accumulated 
   weighted estimate.
   ------------------------------------------ */   
   Fest = 0.00;
   for(j = 0; j < max_idx; j++)
   {
      if(m_pDbase[j].id == id)
	   {  
         Fi = m_pDbase[j].F; 
         Di = 0.00;
         for(i = 0; i < m_NumParams; i++){
            x1 = pX[i];
            x2 = m_pDbase[j].pParams[i];
            Di += ((x2-x1)*(x2-x1));}
	      Di = sqrt(Di);
         Di = 1.00/Di;
         if(Di < Dmin){ 
            Wi = 0.00;}
         else{
            Wi = Di/Dflt;}
         Fest += (Wi*Fi);
      }/* end if() */
   }/* end for() */
      
  return Fest;
}/* end InvDistWSSE() */

/*****************************************************************************
Write()
*****************************************************************************/
void SurrogateDbase::Write(FILE * pFile)
{
}/* end Write() */

/*****************************************************************************
Insert()
*****************************************************************************/
void SurrogateDbase::Insert
(
   ParameterGroup * pGroup, 
   double F, 
   int id, 
   double run_time,
   int mode
)
{
   bool found;
   int min_time_stamp;
   double val, Fmax;
   int i, j, k, max_idx, worst_idx, oldest_idx;   
   static bool bReported = false;
    
   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   if(run_time < 1.00) run_time = 1.00;

   /* -----------------------------------
   Don't store redundant database entries
   Also locate the worst and oldest entries
   for the given model
   ------------------------------------ */ 
   worst_idx = -1;
   oldest_idx = -1;
   for(i = 0; i < max_idx; i++)
   {
      found = true;

      if(m_pDbase[i].id != id)
      {
         found = false;
      }
      else
      {
         //track worst entry
         if(worst_idx == -1){
            worst_idx = i;
            Fmax = m_pDbase[i].F;
         }
         else
         {
            if(m_pDbase[i].F > Fmax)
            {
               worst_idx = i;
               Fmax = m_pDbase[i].F;
            }
         }

         //track oldest entry
         if(oldest_idx == -1){
            oldest_idx = i;
            min_time_stamp = m_pDbase[i].time_stamp;
         }
         else
         {
            if(m_pDbase[i].time_stamp < min_time_stamp)
            {
               oldest_idx = i;
               min_time_stamp = m_pDbase[i].time_stamp;
            }
         }
            
         for(k = 0; k < m_NumParams; k++)
         {
            val = pGroup->GetParamPtr(k)->GetEstVal();
            if(m_pDbase[i].pParams[k] != val)
            {
               found = false;
               break;
            }
         }
      }
      //if a match was found, return without inserting
      if(found == true)
      {
         return;
      }
   }/* end for() */

   //assign insertion index based on mode
   if(m_CurSize < m_MaxSize){
      j = m_CurSize;}
   else if (mode == OVERWRITE_DEFAULT){
      j = (m_CurSize % m_MaxSize);}
   else if(mode == OVERWRITE_OLDEST){
      j = oldest_idx;}
   else if(mode == OVERWRITE_LEAST_FIT){
      j = worst_idx;}

   for(i = 0; i < m_NumParams; i++)
   {
      m_pDbase[j].pParams[i] = pGroup->GetParamPtr(i)->GetEstVal();
   }/* end for() */
   m_pDbase[j].F = F;
   m_pDbase[j].id = id;
   m_pDbase[j].run_time = run_time;
   m_pDbase[j].time_stamp = m_CurSize;
   m_CurSize++;

   if((m_CurSize >= m_MaxSize) && (bReported == false))
   {
      LogError(ERR_ARR_BNDS, "SurrogateDbase::Insert() --> database not large enough to store all model evaluations");
      bReported = true;
   }/* end if() */
}/* end Insert() */

/*****************************************************************************
GetRelativeRunTime()

Compute the run time of the given model id, relative to the maximum computation
time for any of the models.
*****************************************************************************/
double SurrogateDbase::GetRelativeRunTime(int id)
{ 
   int i, j, count, max_idx;
   double avg, max;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   max = 0.00;
   for(i = 0; i < m_NumModels; i++)
   {
      avg = 0.00;
      count = 0;
      for(j = 0; j < max_idx; j++)
      {
         if(m_pDbase[j].id == i) 
         {
            avg += m_pDbase[j].run_time;
            count++;
         }
      }
      m_pAvgRunTimes[i] = avg / count;
      if(m_pAvgRunTimes[i] > max) max = m_pAvgRunTimes[i];
   }

   return m_pAvgRunTimes[id]/max;
}/* end GetRelativeRunTime() */

/*****************************************************************************
GetNumStoredEvals()

Retrieve the number of evals of a given model that are currently stored in the
database.
*****************************************************************************/
int SurrogateDbase::GetNumStoredEvals(int id)
{
   int i, count, max_idx;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   count = 0;
   for(i = 0; i < max_idx; i++)
   {
      if(m_pDbase[i].id == id) 
      {
         count++;
      }
   }

   return count;
}/* end GetNumStoredEvals() */

/*****************************************************************************
LoadBasis()

Load up a radial basis struct with the known values of a given model.

Returns the number of paramter groups loaded into pBasis.
*****************************************************************************/
int SurrogateDbase::LoadBasis(int id, MyPoint * pBasis)
{ 
   int j, count, max_idx;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   count = 0;
   for(j = 0; j < max_idx; j++)
   {
      if(m_pDbase[j].id == id) 
      {
         pBasis[count].F = m_pDbase[j].F;
         pBasis[count].v = m_pDbase[j].pParams;
         pBasis[count].ndim = m_NumParams;
         count++;
      }
   }

   return count;
}/* end LoadBasis() */

/*****************************************************************************
GetBestEntry()

Retrieve the best entry (lowest WSSE) in the database.
*****************************************************************************/
double * SurrogateDbase::GetBestEntry(void)
{
   int j, max_idx;
   double * pBest, F, Fmin;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   pBest = NULL;
   for(j = 0; j < max_idx; j++)
   {
      F = m_pDbase[j].F;
      if(j == 0)
      {
         pBest = m_pDbase[j].pParams;
         Fmin = F;
      }
      else
      {         
         if(F < Fmin)
         {
            Fmin = F;
            pBest = m_pDbase[j].pParams;
         }
      }
   }/* end for() */

   return pBest;
}/* end GetBestEntry() */

/*****************************************************************************
GetBestEntry()

Retrieve the best entry (lowest WSSE) in the database for the given model id.
*****************************************************************************/
DbaseEntry * SurrogateDbase::GetBestEntry(int id)
{
   int j, max_idx;
   DbaseEntry * pBest;
   double F, Fmin;
   bool first = true;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   pBest = NULL;
   for(j = 0; j < max_idx; j++)
   {
      if(m_pDbase[j].id == id)
      {
         F = m_pDbase[j].F;
         if(first == true)
         {
            pBest = &(m_pDbase[j]);
            Fmin = F;
            first = false;
         }
         else
         {         
            if(F < Fmin)
            {
               Fmin = F;
               pBest = &(m_pDbase[j]);
            }
         }
      }/* end if() */
   }/* end for() */

   return pBest;
}/* end GetBestEntry() */

/*****************************************************************************
BcastEntries()

Broadcast database entries to other processors.
*****************************************************************************/
void SurrogateDbase::BcastEntries(ParameterGroup * pGroup, int mode)
{
   int my_id, proc_id, num_procs, j, k, num_params;
   int model_id, max_idx, num_entries;
   double * pData, F, run_time;

   MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   if (num_procs == 1) return;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   num_params = m_NumParams;

   NEW_PRINT("double", num_params+3);
   pData = new double[num_params+3];
   MEM_CHECK(pData);

   for(proc_id = 0; proc_id < num_procs; proc_id++)
   {
      //first broadcast number of entries
      num_entries = max_idx;
      MPI_Bcast(&num_entries, 1, MPI_INTEGER, proc_id, MPI_COMM_WORLD);

      //now broadcast each entry
      for(j = 0; j < num_entries; j++)
      {
         if(my_id == proc_id)
         {
            for(k = 0; k < num_params; k++)
            {
               pData[k] = m_pDbase[j].pParams[k];
            }/* end for() */
            pData[k] = m_pDbase[j].F;
            pData[k+1] = (double)(m_pDbase[j].id);
            pData[k+2] = m_pDbase[j].run_time;
         }/* end if() */
         MPI_Bcast(pData, num_params+3, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
         if(my_id != proc_id)
         {
            pGroup->WriteParams(pData);
            F = pData[num_params];
            model_id = (int)(pData[num_params+1]);
            run_time = pData[num_params+2];
            Insert(pGroup, F, model_id, run_time, mode);
         }/* end if() */
      }/* end for() */
   }/* end for() */

   delete [] pData;
}/* end BcastEntries() */

/*****************************************************************************
BcastBestEntries()

Broadcast the best database entries to other processors.
*****************************************************************************/
void SurrogateDbase::BcastBestEntries(ParameterGroup * pGroup, int mode)
{
   int my_id, proc_id, num_procs, j, k, num_params;
   int model_id, max_idx;
   double * pData, F, run_time;
   DbaseEntry * pBest;

   MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   if (num_procs == 1) return;

   max_idx = m_CurSize;
   if(max_idx > m_MaxSize){max_idx = m_MaxSize;}

   num_params = m_NumParams;

   NEW_PRINT("double", num_params+3);
   pData = new double[num_params+3];
   MEM_CHECK(pData);

   for(proc_id = 0; proc_id < num_procs; proc_id++)
   {
      for(j = 0; j < m_NumModels; j++)
      {
         if(my_id == proc_id)
         {
            pBest = GetBestEntry(j);
//            if(pBest == NULL){
//               printf("couldn't get best entry for %d\n", j);}
//            printf("proc %d, pBest[0] = %lf\n", proc_id, pBest->pParams[0]);
            for(k = 0; k < num_params; k++)
            {
               pData[k] = pBest->pParams[k];
            }/* end for() */
            pData[k] = pBest->F;
            pData[k+1] = (double)(pBest->id);
            pData[k+2] = pBest->run_time;
         }/* end if() */
//         printf("proc %d : broadcasting best model %d, pData = %0x\n", proc_id, j, (unsigned int)pData);
         MPI_Bcast(pData, num_params+3, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
//         printf("proc %d : done broadcasting best model %d, pData = %0x\n", proc_id, j, (unsigned int)pData);
         if(my_id != proc_id)
         {
            pGroup->WriteParams(pData);
            F = pData[num_params];
            model_id = (int)(pData[num_params+1]);
            run_time = pData[num_params+2];
            Insert(pGroup, F, model_id, run_time, mode);
         }/* end if() */
      }/* end for() */
   }/* end for() */

   delete [] pData;
}/* end BcastBestEntries() */
