/******************************************************************************
File      : LatinHypercube.cpp
Author    : L. Shawn Matott
Copyright : 2006, L. Shawn Matott

Encapsulates a lating hypercube sampling strategy for initializing populations.

Version History
06-14-06    lsm   created
******************************************************************************/
#include <stdio.h>
#include <math.h>

#include "LatinHypercube.h"

#include "Exception.h"
#include "Utility.h"
#include "StatUtility.h"

/******************************************************************************
Destroy()
******************************************************************************/
void LatinHypercube::Destroy(void)
{
   int i;

   delete [] m_pCount;
   for(i = 0; i < m_Rows; i++){ delete [] m_pVals[i];}
   delete [] m_pVals;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR
******************************************************************************/
LatinHypercube::LatinHypercube(int rows, int cols)
{
   int i, j;

   NEW_PRINT("double *", rows);
   m_pVals = new double *[rows];
   MEM_CHECK(m_pVals);

   NEW_PRINT("int", rows);
   m_pCount = new int [rows];
   MEM_CHECK(m_pCount);

   for(i = 0; i < rows; i++)
   {
      NEW_PRINT("double", cols);
      m_pVals[i] = new double[cols];
      MEM_CHECK(m_pVals[i]);

      m_pCount[i] = cols;
   }/* end for() */

   for(i = 0; i < rows; i++)
   {
      for(j = 0; j < cols; j++)
      {
         m_pVals[i][j]  = 0.00;
      }
   }

   m_Rows = rows;
   m_MaxCols = m_Cols = cols;

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
ReDim()

Redimension the number of columns. Cannot increase, but can decrease.
******************************************************************************/
void LatinHypercube::ReDim(int cols)
{
   if(cols > m_MaxCols)
   {
      LogError(ERR_ARR_BNDS, "Can't redimension hypercube");
      ExitProgram(1);
   }
   m_Cols = cols;
} /* end default ReDim() */

/******************************************************************************
InitRow()

Initialize a row of the hypercube sampling matrix using a uniform distribtion.
******************************************************************************/
void LatinHypercube::InitRow(int row, double min, double max)
{
   int i;
   double step, lwr;

   step = (max - min)/m_Cols;
   lwr = min;
   for(i = 0; i < m_Cols; i++)
   {
      m_pVals[row][i] = lwr + step*((double)MyRand()/(double)MY_RAND_MAX);      
      lwr += step;
   }
   m_pCount[row] = m_Cols; //reset sample count
} /* end uniform InitRow() */

/******************************************************************************
InitRow()

Initialize a row of the hypercube sampling matrix using a Gaussian (i.e Normal) 
distribution. The distribution is truncated so that samples lie between min 
and max. The truncated distribution is split into m_Cols intervals of equal
probability. Each interval is then randomly sampled.
******************************************************************************/
void LatinHypercube::InitRow(int row, double min, double max, double sd)
{
   int i;
   double z_min, z_max; //std. normal min and max
   double p_min, p_max; //cumulative probabilities
   double p_step; //step size in probability units
   double z_step; //step size in std. normal units
   double z_rand; 
   double avg;

   avg = 0.5*(max+min);
   z_max = (max - avg)/sd;
   z_min = (min - avg)/sd;
   p_max = StdNormCDF(z_max);
   p_min = StdNormCDF(z_min);
   p_step = (p_max - p_min)/ m_Cols;
   
   for(i = 0; i < m_Cols; i++)
   {
      z_max = StdNormInvCDF(p_min + p_step);
      z_step = z_max - z_min;
      z_rand = z_min + z_step*((double)MyRand()/(double)MY_RAND_MAX);

      m_pVals[row][i] = avg + sd*z_rand;
      p_min += p_step;
      z_min = z_max;
   }
   m_pCount[row] = m_Cols; //reset sample count
} /* end default InitRow() */

/******************************************************************************
SampleRow()

Sample from the hypercube matrix.
******************************************************************************/
double LatinHypercube::SampleRow(int row)
{
   double sample;
   int i;

   i = (MyRand() % m_pCount[row]);
   sample = m_pVals[row][i];
   
   //reorder list
   for(i = i; i < (m_Cols-1); i++)
   {
      m_pVals[row][i]  = m_pVals[row][i+1];
   }
   m_pVals[row][m_Cols-1] = sample;
   m_pCount[row]--;

   return sample;
} /* end Sample() */

