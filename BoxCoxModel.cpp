/******************************************************************************
File     : BoxCoxModel.cpp
Author   : L. Shawn Matott
Copyright: 2013, L. Shawn Matott

Computes an objective function value for a given BoxCox transformation. The obj.
function value measures the degree to which the transformation incurs normality
on the transformed and weighted residuals.

Version History
01-09-13    lsm   Created
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "BoxCoxModel.h"

#include "StatUtility.h"
#include "Utility.h"

double BoxCoxVal(double y, double w, double c);
double BoxCoxNPP(double * v, int n);

int BoxCoxModel(void)
{  
   FILE * pIn = fopen(BOX_IN_FILE, "r");
   FILE * pOut = fopen(BOX_OUT_FILE, "w");
   if((pIn == NULL) && (pOut == NULL)) return -1;
   if(pIn == NULL)
   {
      fprintf(pOut, "Unable to open input file |%s|\n", BOX_IN_FILE);
      fclose(pOut);
      return -1;
   }
   if(pOut == NULL)
   {
      printf("Unable to open output file |%s|\n", BOX_OUT_FILE);
      fclose(pIn);
      return -1;
   }

   //simple file parsing, easily broken if input file is missing items or contains typos
   double p = 1.00;
   int n = 0;
   fscanf(pIn, "LAMBDA=%lf\n", &p);
   fscanf(pIn, "NUM_DATA_POINTS=%d\n", &n);
   double x, y, w, *r;
   r = new double[n];
   for(int i = 0; i < n; i++)
   {
      x=y=w=1.00;
      fscanf(pIn, "%lf %lf %lf\n", &x, &y, &w);
      r[i]=BoxCoxVal(x,w,p)-BoxCoxVal(y,w,p);
   }/* end for() */
   fclose(pIn);

   double npp = BoxCoxNPP(r, n);

   fprintf(pOut,"ObjFunc=%E\n",-npp);
   fprintf(pOut,"R-squared=%lf\n",npp*npp);
   fprintf(pOut,"NPP=%lf\n",npp);
   fprintf(pOut, "LAMBDA=%lf\n", p);
   fprintf(pOut, "NUM_DATA_POINTS=%d\n", n);
   fprintf(pOut,"RESIDUALS\n");
   for(int i = 0; i < n; i++)
      fprintf(pOut, "%d\t%E\n", i, r[i]);

   fclose(pOut);
   delete [] r;
   return (0);
} /* end BoxCoxModel() */


/******************************************************************************
BoxCoxVal()

Perform Box-Cox transformation on the given input value.
******************************************************************************/
double BoxCoxVal(double y, double w, double c)
{  
  y=w*y; //apply weight

  //y must be positive, if not: don't perform transformation and log error
  if(y < 0.00)
  {
    printf("Couldn't perform Box-Cox transformation, data is non-positive!");
    return y;
  }

  double h = 0.00;
  if(c != 0)
  {
    h = (pow(y, c) - 1.00)/c;
  }
  else //natural log transformation
  {
    h = log(y);
  }
  return h;
} /* end BoxCoxVal() */

/******************************************************************************
BoxCoxNPP()

Calculates the normal probability plot correlation coefficient for a set of
values.
******************************************************************************/
double BoxCoxNPP(double * v, int n)
{      
   int i;
   double pi;
   double * w = new double[n]; //expected values

   SortInc(v, n);
   double mean = CalcMean(v, n);

   /*----------------------------------------------
   Compute the expected values of std. norm. order 
   stats. using the Snedecor & Cochran approximation
   described by David W. Sabo (BCIT) in "Normal 
   Probability Plots", page #3. 
   ----------------------------------------------*/
   for(i = 0; i < n; i++)
   {
      pi = (double(i+1) - 0.375)/((double)n + 0.25);
      w[i] = StdNormInvCDF(pi);
   }/* end for() */   
   
   /*----------------------------------------------
   Compute the numerator of Equation 5 of "Methods 
   and Guidelines for Effective Model Calibration",
   page 23.
   ----------------------------------------------*/
   double numer = 0.00;
   for(i = 0; i < n; i++)
      numer += ((v[i] - mean)*w[i]);
   numer *= numer;

   /*----------------------------------------------
   Calc. the denominator of Equation 5 of "Methods 
   and Guidelines for Effective Model Calibration",
   page 23.
   ----------------------------------------------*/
   double denom1 = 0.00;
   double denom2 = 0.00;
   for(i = 0; i < n; i++)
   {
      denom1 += ((v[i] - mean)*(v[i] - mean));
      denom2 += (w[i]*w[i]);
   }/* end for() */   
   
   double npp = (numer / (denom1 * denom2));

   delete [] w;

   return sqrt(npp);
}/* end BoxCoxNPP() */

