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
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>

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

   for(i = 0; i < m_Rows; i++){ delete [] m_pVals[i];}
   delete [] m_pVals;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR
******************************************************************************/
LatinHypercube::LatinHypercube(int numberOfSamples, int numberOfParameters, int numberOfBins)
{
   int i, j;

   // Set the values to hold for future analyses
   m_Rows = numberOfSamples;
   m_Cols = numberOfParameters;
   m_MaxCols = numberOfBins;

   // Adjust the to make sure the sample has sufficient combinations when sampleing
   while (m_MaxCols < m_Rows) {
       m_MaxCols *= 2;
   }

   // Create the sample initial array
   m_pVals = new double *[m_Rows];
   MEM_CHECK(m_pVals);

   // Create the sample subarrays
   for(i = 0; i < m_Rows; i++) {
      // Initialize the array
      m_pVals[i] = new double[m_Cols];
      MEM_CHECK(m_pVals[i]);

      // Set the value to zero
      for (j = 0; j < m_Cols; j++) {
          m_pVals[i][j] = 0.00;
      }
   }

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
CreateUniformSample()

Initializes the hypercube sampling matrix using a uniform distribtion.
******************************************************************************/
void LatinHypercube::CreateUniformSample(double* lowerBounds, double *upperBounds) {

    // Divide the sample space into bins
    double* binWidth = new double[m_Cols];

    for (int i = 0; i < m_Cols; i++) {
        binWidth[i] = (upperBounds[i] - lowerBounds[i]) / (double)(m_MaxCols-1);
    }

    // Create a vector containing the bin numbers for all the inputs for sampling.
    // This is more efficient that randomly sampling and checking against the counter
    // Create the parent vector
    std::vector<std::vector<double>> availableIndices;

    // Fill the vector
    for (int i = 0; i < m_Cols; i++) {
        // Create the temporary vector
        std::vector<double> temp;

        // Fill it with values
        for (int j = 0; j < m_MaxCols; j++) {
            temp.push_back((double)j);
        }

        // Obtain a time-based seed
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        // Shuffle the indices to simplify later indexing
        std::shuffle(temp.begin(), temp.end(), std::default_random_engine(seed));

        // Push the temporary vector to the parent vector
        availableIndices.push_back(temp);
    }

    // Begin selecting the bins and filling the rows
    for (int i = 0; i < m_Rows; i++) {
        for (int j = 0; j < m_Cols; j++) {
            // Set the value into the function
            m_pVals[i][j] = availableIndices[j].back() * binWidth[j] + lowerBounds[j];

            // Remove the last index of the bin
            availableIndices[j].pop_back();
        }
    }
} 

/******************************************************************************
InitRow()

Initialize a row of the hypercube sampling matrix using a Gaussian (i.e Normal) 
distribution. The distribution is truncated so that samples lie between min 
and max. The truncated distribution is split into m_Cols intervals of equal
probability. Each interval is then randomly sampled.
******************************************************************************/
//void LatinHypercube::InitRow(int row, double min, double max, double sd)
//{
//   int i;
//   double z_min, z_max; //std. normal min and max
//   double p_min, p_max; //cumulative probabilities
//   double p_step; //step size in probability units
//   double z_step; //step size in std. normal units
//   double z_rand; 
//   double avg;
//
//   avg = 0.5*(max+min);
//   z_max = (max - avg)/sd;
//   z_min = (min - avg)/sd;
//   p_max = StdNormCDF(z_max);
//   p_min = StdNormCDF(z_min);
//   p_step = (p_max - p_min)/ m_Cols;
//   
//   for(i = 0; i < m_Cols; i++)
//   {
//      z_max = StdNormInvCDF(p_min + p_step);
//      z_step = z_max - z_min;
//      z_rand = z_min + z_step*((double)MyRand()/(double)MY_RAND_MAX);
//
//      m_pVals[row][i] = avg + sd*z_rand;
//      p_min += p_step;
//      z_min = z_max;
//   }
//   m_pCount[row] = m_Cols; //reset sample count
//} /* end default InitRow() */


