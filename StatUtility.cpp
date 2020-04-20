/******************************************************************************
File     : StatUtility.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

StatUtility contains C routines that assist in statistical calculations.

Version History
06-19-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   replaced magic # with NEARLY_ZERO
01-01-07    lsm   Added Durbin-Watson and runs tests.
******************************************************************************/
#include <string.h>
#include <math.h>

#include "StatUtility.h"
#include "Utility.h"
#include "MyErf.h"
#include "Exception.h"

/******************************************************************************
GammaLn()

Compute natural log of the gamma function.
******************************************************************************/
double GammaLn(double val)
{
   int j;
   double x, y, tmp, ser;
   
   static double coeff[6] = 
   { 76.18009172947146, -86.50532032941677, 24.01409824083091, 
    -1.231739572450155, 0.001208650973866179, -0.5395239384953e-5};

   y = x = val;
   tmp = x + 5.5;
   tmp -= (x + 0.5) * log(tmp);
   ser = 1.000000000190015;

   for(j = 0; j < 6; j++)
   {
      y++;
      ser += coeff[j] / y;
   }/* end for() */       

   return (-tmp + log(2.5066282746310005 * ser / x));
}/* end GammaLn() */

/******************************************************************************
CalcStdDev()

Returns the standard deviation of a list of numbers.
******************************************************************************/
double CalcStdDev(Ironclad1DArray v, int size, int ctType)
{
   int i;
   double mean;
   double sum;
   double v1;
   double tmp;
   int ncensored = 0;

   if(ctType == CENTRAL_TEND_PCTILE)
   {
      // estimate std dev using the raw data
      // half the difference between 84th percentile and 16th percentile
      double * vv = new double[size];
      memcpy(vv, v, size*sizeof(double));
      mean = CalcMedian(vv, size);
      tmp = vv[(int)(size*0.84)];
      tmp -= vv[(int)(size*0.16)];
      tmp *= 0.5;
      delete [] vv;
      return tmp;
   }
   else if(ctType == CENTRAL_TEND_MEAN)
   {
      mean = CalcMean(v, size);
   }
   else //median
   {
      double * vv = new double[size];
      memcpy(vv, v, size*sizeof(double));
      mean = CalcMedian(vv, size);
      delete [] vv;
   }

   sum = 0.00;

   for(i = 0; i < size; i++)
   {
      v1 = v[i];
      tmp = ((v1 - mean) * (v1 - mean));
      if(CheckOverflow(tmp) == true)
         ncensored++;
      else
         sum += tmp;
   } /* end for() */
  
   return sqrt(sum / (double)(size - ncensored));
} /* end CalcStdDev() */

/******************************************************************************
CalcSkewness()

Returns the skewness of a list of numbers.
******************************************************************************/
double CalcSkewness(Ironclad1DArray v, int size)
{
   int i;
   double n = (double)size;
   double mean, sd;
   double sum;
   double v1;

   mean = CalcMean(v, size);
   sd = CalcStdDev(v, size,CENTRAL_TEND_MEAN);

   sum = 0.00;

   for(i = 0; i < size; i++)
   {
      v1 = (v[i] - mean)/sd;
      sum += (v1*v1*v1);
   } /* end for() */
  
   return (n/((n-1)*(n-2)))*sum;
}/* end CalcSkewness() */

/******************************************************************************
CalcKurtosis()

Returns the kurtosis of a list of numbers.
******************************************************************************/
double CalcKurtosis(Ironclad1DArray v, int size)
{
   int i;
   double n = (double)size;
   double mean, sd;
   double sum;
   double v1;

   mean = CalcMean(v, size);
   sd = CalcStdDev(v, size,CENTRAL_TEND_MEAN);

   sum = 0.00;

   for(i = 0; i < size; i++)
   {
      v1 = (v[i] - mean)/sd;
      sum += (v1*v1*v1*v1);
   } /* end for() */
  
   return (((n*(n+1))/((n-1)*(n-2)*(n-3)))*sum - (3*(n-1)*(n-1))/((n-2)*(n-3)));
}/* end CalcKurtosis() */

/******************************************************************************
GetCritValNormPPCC()

Returns the 0.05 significance level critical value for the Normal Probability
Plot Correlation Coefficient. The value is interpolated from tables provided 
here: http://www.itl.nist.gov/div898/handbook/eda/section3/eda3676.htm

The input parameter n is the sample size, which must be positive.
******************************************************************************/
double GetCritValNormPPCC(int n)
{
  if(n < 0) n = 0;
  double c;
  //best-fit interpolation constants for a dual-langmuir type of fit
  double G2 = 0.216924383;
  double G3 = 0.164415744;
  double G4 = -2.21E-04;
  double G5 = -0.327616476;
  double G6 = 1.42119E-06;
  double G7 = 0.781698141;
  double G8 = 1.03399226;

  c=(G2*pow((G3*n),G8))/(1+pow((G3*n),G8))+(G4*(G5*n))/(1+(G5*n))+G6*n+G7;

  if(c > 1.00) c = 1.00;
  return c;
}/* end GetCritValNormPPCC() */

/******************************************************************************
CalcMean()

Returns the mean of a list of numbers
******************************************************************************/
double CalcMean(Ironclad1DArray v, int size)
{
  int i;
  double sum;

  sum = 0.00;
  for(i = 0; i < size; i++){ sum += v[i];}
  return (sum / (double)size);
} /* end CalcMean() */

/******************************************************************************
AutorunFunctionTest()

Compute the autorun function and assocaited critical values. It is a measure of
autocorrelation described in:
   McKenzie E. 1984. The Autorun Function: A Non-Parameteric Autocorrelation 
   Function. Journal of Hydrology, vol. 67, no. 1-4, pg. 45-53.
******************************************************************************/
void AutorunFunctionTest(double * residuals, int num, 
                         double * r1,  // lag-1 autorun function 
                         double * var, // variance of r1 
                         double * vpx, // approximate variance of r1 
                         double * med, // median of residuals
                         int * nSur, // number of surplues (values > median)
                         int * nDef, // number of deficits (values <= median)
                         int * n1, // number of lag-1 surpluses
                         double * clwr, // lower-tail critical value
                         double * cupr) // upper-trail critical values
{
   int i;
   double * pR = new double[num]; //scratch vector for residuals

   //determine median, noting that call will sort values
   for(i = 0; i < num; i++) pR[i] = residuals[i];
   *med = CalcMedian(pR, num);

   //count surpluses and deficits
   *nSur = *nDef = 0;
   for(i = 0; i < num; i++)
   {  
      if(pR[i] > (*med)) (*nSur)++;
      else (*nDef)++;
   }

   //count lag-1 surpluses
   *n1 = 0;
   for(i = 1; i < num; i++)
   {  
      if((residuals[i] > (*med)) && (residuals[i-1] > (*med))) (*n1)++;
   }

   double m = *med;
   double N = (double)num;
   double nk = (double)(*n1);
   double Ek = (double)(*n1);
   double  k = 1.00;
   double a, b, c, sd;

   //from page 46 of McKenzie(1984) --- last sentance of 1st full paragraph 
   *r1 = (2.00*nk)/(N - k);

   //from page 47 of McKenzie(1984)
   Ek = (N - 2)/(2*(N - 1));

   //from equation (3) of McKenzie(1984)
   a = -(2*N*N - 9*N + 6);
   b = (N*(N - 1)*(7*N - 26))*0.5;
   c = -N*N*(N - 1)*(N - 4);

   //from equation (3) of McKenzie(1984), but with corrected denominator
   *var = (a + (b/(N - k)) + (c/((N - k)*(N- k))))/(2*(N - 1)*(N - 1)*(N - 3));

   //from page 47 of McKenzie(1984)
   *vpx = 1/(4*(N - 3*k));

   sd = sqrt(*var);

   //critical values at the alpha=0.1 significance level
   *clwr = Ek - 1.645*sd;
   *cupr = Ek + 1.645*sd;

   delete [] pR;
}/* end AutorunFunctionTest() */

/******************************************************************************
CalcMedian()

Returns the median of a list of numbers
******************************************************************************/
double CalcMedian(Unmoveable1DArray v, int size)
{
  int i = size/2;

  SortInc(v,size);

  //even
  if((size % 2) == 0){ return 0.5*(v[i]+v[i-1]); }
  //odd
  return v[i];
} /* end CalcMedian() */

/******************************************************************************
FdistPDF()

Calculates the probability density function of the F-distribution. 

Code is based on formula presented in "Applied Statistics and Probability 
for Engineers", Montgomery and Runger, 1994. Page 316, Theorem 6-5.
******************************************************************************/
double FdistPDF(int u, int v, double x)
{
   double u1 = (double)u;
   double v1 = (double)v;
   double val;   

   double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

   if(x < 0) { return 0;}   
   
   tmp1 = GammaLn((u1 + v1) / 2.00);
   tmp2 = GammaLn(u1 / 2.00);
   tmp3 = GammaLn(v1 / 2.00);
   tmp4 = exp(tmp1 - (tmp2 + tmp3));
   tmp5 = pow((u1 / v1), (u1 / 2.00));
   tmp6 = tmp4 * tmp5 * pow(x, ((u1 / 2.00) - 1.00));
   tmp7 = ((u1 / v1) * x) + 1.00;
   tmp7 = pow(tmp7, ((u1 + v1) / 2.00));

   val = (tmp6 / tmp7);

   return val;
}/* end FdistPDF() */

/******************************************************************************
FdistCDF()

Calculates the cumulative density function of the F-distribution. The args are
u = d.o.f. in numerator
v = d.o.f. in denominator
xLwr = lower limit
xUpr = upper limit

Returns the probability that a F-distributed variable with u and v degress of 
freedom is less than or equal to xUpr and greter than or equal to xLwr 
(i.e. P[xLwr <= X <= xUpr]).
******************************************************************************/
double FdistCDF(int u, int v, double xLwr, double xUpr)
{
   int j;
   double lwr, upr, mid; //integration limits
   double dt;
   double sum, old;
   double Flwr, Fupr, Fmid;
   double stop_val = 1E-6;
   
   if(v <= 0)
   {   
      LogError(ERR_BAD_DOF, "Not enough degrees of freedom for FdistCDF()");
      return -1.00;
   }

   /*--------------------------------------
   determine integration limits.
   ---------------------------------------*/
   if(xLwr <= 0.00){xLwr = 0.00;}
   if(xUpr <= 0.00){return 0.00;}
   if(xLwr == xUpr){return 0.00;}
   if((xLwr == 0.00) && (xUpr == 1.00) && (u == v)){return 0.50;}

   if(u == 1)
   { 
      Fupr = (2.00*StudentCDF(v, sqrt(xUpr)) - 1.00);
      Flwr = (2.00*StudentCDF(v, sqrt(xLwr)) - 1.00);

      return (Fupr - Flwr);
   }/* end if() */
   
   lwr = xLwr; 
   upr = xUpr;

   /*-----------------------------------------------------
   integrate over the limits, the integration technique
   is an iterative form of the trapezoidal rule. The 
   number of intervals is doubled each iteration until 
   the area calculation has converged.
   ------------------------------------------------------*/
   Flwr = FdistPDF(u, v, lwr);
   Fupr = FdistPDF(u, v, upr);
   dt = (upr - lwr);
   //initially apply trapezoidal rule at the end points
   sum = 0.5*dt*(Flwr+Fupr); 
   old = sum;
   j = 0;   
   //iterate until converged or at least 5 times
   while((fabs(sum - old) > stop_val) || (j < 5))
   {
      j++;      
      dt /= 2; //halve step size
      mid = lwr + dt; //initial value of points
   
      old = sum;
      sum = 0.00;
      //calculate F() of points
      while(mid <= upr)
      {
         Fmid = FdistPDF(u, v, mid);
         sum +=  Fmid;
         //a step of 2*dt skips previously computed points    
         mid += (2.00*dt); 
      }/* end while() */
      sum *= dt;      
      sum += (old * 0.5);
   }/* end while() */   
      
   return (sum);
}/* end FdistCDF() */

/******************************************************************************
FdistInvCDF()

Calculates the F-distribution upper-tail percentage point, (i.e. given  a 
probability p, compute the value of x such that the P[X < x] = p).
******************************************************************************/
double FdistInvCDF(int u, int v, double p)
{
   int j;
   double x, upr, lwr;
   double F;
   double stopVal = 1E-6;

   if(v <= 0)
   {   
      LogError(ERR_BAD_DOF, "Not enough degrees of freedom for FdistInvCDF()");
      return -1.00;
   }

   if(p <= 0.00){ return 0.00;}
   if(p >= 1.00){ p = 1.00;} //max. probability is 1

   //initial guess is at center of distribution
   x = 1.00;
   F = FdistCDF(u, v, 0.00, x);

   upr = 1.00;
   lwr = 0.00;
   
   //find upper limit
   while(F < p) 
   { 
      lwr = x;
      x *= 2;
      //add in the area between lwr and x
      F += FdistCDF(u, v, lwr, x); 
   }/* end while() */   

   //converge on answer
   j = 0;
   while((fabs(F - p) > stopVal) || (j < 5))
   {
      if(F >= p) 
      { 
         upr = x;
         x = 0.5*(upr + lwr);
         //subtract the area between x and upr
         F -= FdistCDF(u, v, x, upr);
      }/* end if() */
      else 
      { 
         lwr = x;
         x = 0.5*(upr + lwr);
         //add the area between lwr and x
         F += FdistCDF(u, v, lwr, x);
      }/* end else() */

      j++;
   }/* end while() */
   
   return x;
}/* end FdistInvCDF() */

/******************************************************************************
StudentPDF()

Calculates the probability density function of the Student's t-distribution, 
given dof (degrees of freedom) and x.
******************************************************************************/
double StudentPDF(int dof, double x)
{
   double coeff; //coefficient, only depends on dof
   double e; //exponent, only depends on dof
   double tmp1, tmp2, g1, g2;
   double val;

   /*-----------------------------------
   Calculate terms that rely only on dof
   -----------------------------------*/
   tmp1 = 0.5*((double)dof + 1.00);
   e = -tmp1; 
   tmp2 = 0.5*(double)dof;
   g1 = GammaLn(tmp1);
   g2 = GammaLn(tmp2);
   tmp1 = exp(g1 - g2);
   coeff = tmp1 / sqrt((double)dof*MY_PI);

   /*-----------------------------------
   calculate the x-dependent term
   ------------------------------------*/
   tmp1  = (1.00 + ((x*x)/(double)dof));
   tmp2 = pow(tmp1, e);
   val = (coeff*tmp2);
      
   return (val);
}/* end StudentPDF() */

/******************************************************************************
StudentCDF()

Calculates the Student's cumulative t-distribution, given dof (degrees of 
freedom) and x.

Returns the probability that a t-distributed variable with dof degress of 
freedom is less than or equal to x (i.e. P[X <= x]).
******************************************************************************/
double StudentCDF(int dof, double x)
{
   double P; //desired probability
   double lwr, upr, mid; //integration limits
   double dt;
   double sum, old;
   double Flwr, Fupr, Fmid;   
   double stop_val = 1E-6;
   int j;

   //area under PDF from -infinity to 0 is 0.5
   P = 0.5; 

   /*--------------------------------------
   determine integration limits.
   ---------------------------------------*/
   if(x == 0.00) { return P;}
   else if(x > 0.00){lwr = 0.00; upr = x;}
   else {lwr = x; upr = 0.00;}

   /*-----------------------------------------------------
   integrate over the limits, the integration technique
   is an iterative form of the trapezoidal rule. The 
   number of intervals is doubled each iteration until 
   the area calculation has converged.
   ------------------------------------------------------*/
   Flwr = StudentPDF(dof, lwr);
   Fupr = StudentPDF(dof, upr);
   dt = (upr - lwr);
   //initially apply trapezoidal rule at the end points
   sum = 0.5*dt*(Flwr+Fupr); 
   old = sum;
   j = 0;   
   //iterate until converged or at least 5 times
   while((fabs(sum - old) > stop_val) || (j < 5))
   {
      j++;
      dt /= 2; //halve step size
      mid = lwr + dt; //initial value of points
   
      old = sum;   
      sum = 0.00;
      //calculate F() of points
      while(mid <= upr)
      {
         Fmid = StudentPDF(dof, mid);
         sum +=  Fmid;
         //a step of 2*dt skips previously computed points    
         mid += (2.00*dt); 
      }/* end while() */
      sum *= dt;
      sum += (old * 0.5);
   }/* end while() */

   if(x > 0.00) { P += sum;}
   else { P -= sum;}
      
   return (P);
}/* end StudentCDF() */

/******************************************************************************
StudentInvCDF()

Calculates the student-distribution upper-tail percentage point, (i.e. given  a 
probability p, compute the value of x such that P[X < x] = p).
******************************************************************************/
double StudentInvCDF(int dof, double p)
{
   double x, upr, lwr;
   double F;
   double stopVal = 1E-6;
   bool flip = false;

   if(p <= 0.00){ p = 0.00;} //min. probability is 0
   if(p >= 1.00){ p = 1.00;} //max. probability is 1
   if(p == 0.50){ return 0.00;}   
   if(p < 0.50){ flip = true; p = 1 - p;}

   //initial guess is at center of distribution
   x = 0.00; F = 0.50;

   upr = 1.00; lwr = 0.00;
   
   //find upper limit
   while(F < p) 
   { 
      lwr = x;
      x = upr;
      upr *= 2;
      F = StudentCDF(dof, x);
   }/* end while() */   

   //converge on answer
   while(fabs(F - p) > stopVal)
   {      
      if(F >= p) { upr = x;}   
      else { lwr = x;}

      x = 0.5*(upr + lwr);
      F = StudentCDF(dof, x);
   }/* end while() */
   
   if(flip == true){ return -x;}
   return x;
}/* end StudentInvCDF() */

/******************************************************************************
StdNormPDF()

Calculates the probability density function of the standard normal distribution.
This is a normal distribution with mean = 0 and std. dev = 1.
******************************************************************************/
double StdNormPDF(double x)
{
   double val;
   
   val = (exp(-(x*x)*0.5)/(sqrt(2.00 * MY_PI)));
   
   return (val);
}/* end StdNormPDF() */

/******************************************************************************
StdNormCDF()

Calculates the cumulative density function of the standard normal distribution.
This is a normal distribution with mean = 0 and std. dev = 1.

Returns the probability that a normally variable is less than or equal to x 
(i.e. P[X <= x]).

NOTE: Utilizes the relationship of the std. normal distribution to the error
function.
******************************************************************************/
double StdNormCDF(double x)
{
   double val;
   
   val = 0.5*(MyErf(x/sqrt(2.00)) + 1.00);
   
   return (val);
}/* end StdNormCDF() */

/******************************************************************************
StdNormInvCDF()

Calculates the standard normal distribution upper-tail percentage point, 
(i.e. given  a probability p, compute the value of x such that P[X < x] = p).
******************************************************************************/
double StdNormInvCDF(double p)
{
   double x, upr, lwr;
   double F;
   double stopVal = NEARLY_ZERO;
   bool flip = false;

   if(p <= 0.00){ p = 0.00;} //min. probability is 0
   if(p >= 1.00){ p = 1.00;} //max. probability is 1
   if(p == 0.50){ return 0.00;}
   if(p < 0.50){ flip = true; p = 1 - p;}

   //initial guess is at center of distribution
   x = 0.00; F = 0.50;
   lwr = 0.00; upr = 1.00;
   
   //find upper limit
   while(F < p) 
   { 
      lwr = x;
      x = upr;
      upr *= 2;

      F = StdNormCDF(x);
      if(fabs(F - p) < stopVal){ break;}
   }/* end while() */   

   //converge on answer
   while(fabs(F - p) > stopVal)
   {      
      if(F > p) { upr = x;}   
      else { lwr = x;}

      x = 0.5*(upr + lwr);
      F = StdNormCDF(x);
   }/* end while() */

   if(flip == true){ return -x;}   
   return x;
}/* end StdNormInvCDF() */

/******************************************************************************
RunsTest()

Perform a runs test for autocorrelation of residuals. The residuals should be 
placed in an appropriate order by the calling routine (e.g. ordered by
increasing value of observation).

Code is based on the formulation presented by Swed and Eisenhart in "Tables for 
testing randomness of grouping in a sequence of alternatives", The Annals of 
Mathematical Statistics, vol. 14, no. 1, March 1943, pg. 66-87.

NOTE: a value of zero is counted as positive.

Input:
   residuals : Array of residuals, runs test examines the sign of the residuals.
   n : size of residuals array
Output:
   Returns true if successful (npos and nneg > 1) and sets the following output 
   values:
      nPos   : number of positive values
      nNeg   : number of negative values
      nRuns  : number of runs
      clwr   : critical value for lower tail (nRuns < clwr)
      cupr   : critical value for upper tail (nRuns > cupr)
   Return false if npos or nneg is less than or equal to 1 (too few samples)
******************************************************************************/
bool RunsTest(double * residuals, int num, int * nPos, int * nNeg, int * nRuns, 
              int * clwr, int * cupr)
{
   int i, s1, s2;
   
   *nPos = *nNeg = 0; 
   *nRuns = 1; 

   for(i = 0; i < num; i++)
   {
      if(residuals[i] < 0.00) (*nNeg)++;
      else                    (*nPos)++;
   }
   if((*nPos < 2) || (*nNeg < 2)){*clwr = *cupr = 0; return false;}

   for(i = 1; i < num; i++)
   {
      if(residuals[i] < 0.00) s1 = -1;
      else s1 = +1;
      if(residuals[i-1] < 0.00) s2 = -1;
      else s2 = +1;
      if(s1 != s2) (*nRuns)++;
   }/* end for() */

   /* -------------------------------------------------------------------------
   Create a Runs Test CDF table ---> for deubgging purposes...
   ------------------------------------------------------------------------- */
   /*
   FILE * pRT = fopen("OstRunsTestTable.txt", "w");
   fprintf(pRT, "alpha  nruns\n");
   for(i=0; i<(*nPos + *nNeg); i++)
     fprintf(pRT, "%5lf  %d\n", RunsTestCDF(i, *nPos, *nNeg), i);
   fclose(pRT);
   */
   /* -----------------------------------------------------------------------
   Spot check results against tabulated values (see Sewd and Eisenhart)
   Note ---- for debugging purposes only.
   ----------------------------------------------------------------------- */
   /*
   int t1 = InvRunsTestCDF(0.025, 4, 9);
   int t2 = InvRunsTestCDF(0.975, 4, 9);

   t1 = InvRunsTestCDF(0.025, 20, 9);
   t2 = InvRunsTestCDF(0.975, 20, 9);

   t1 = InvRunsTestCDF(0.01, 13, 16);
   t2 = InvRunsTestCDF(0.99, 13, 16);

   t1 = InvRunsTestCDF(0.005, 18, 19);
   t2 = InvRunsTestCDF(0.995, 18, 19);

   t1 = InvRunsTestCDF(0.05, 38, 38);
   t2 = InvRunsTestCDF(0.95, 38, 38);

   double t3 = RunsTestCDF(4, 2, 15);
   t3 = RunsTestCDF(4, 4, 11);
   t3 = RunsTestCDF(5, 3, 17);
   */
     
   //determine critical values
   *clwr = InvRunsTestCDF(0.05, *nPos, *nNeg);
   *cupr = InvRunsTestCDF(0.95, *nPos, *nNeg);

   return true;
}/* end RunsTest() */

/******************************************************************************
InvRunsTestCDF()

Determine the number of runs required to satisfy the given significance level.
******************************************************************************/
int InvRunsTestCDF(double a, int m, int n)
{
   int i;
   for(i = 0; i < (m+n); i++)
   {
     if(a <= 0.5)
     {
       if(RunsTestCDF(i, m, n) > a)
       {
          return (i-1);
       }
     }
     else
     {
       if(RunsTestCDF(i, m, n) > a)
       {
          return (i);
       }
     }
   }/* end for() */
   return (i-1);
}/* end InvRunsTestCDF() */

/******************************************************************************
RunsTestCDF()

Computes the p-value for a given number of runs if there are m members in 
"Group A" and n members in "Group B".
******************************************************************************/
double RunsTestCDF(int nRuns, int m, int n)
{
   double p, fu, sum;
   int k, u, tmp;
   //ensure m <= n
   if(m > n) { tmp = m; m = n; n = tmp; }

   sum = 0.00;
   for(u = 2; u <= nRuns; u++)
   {
      //u is even
      if((u % 2) == 0)
      {
         k = u / 2;
         fu = (2 * nCr(m-1, k-1) * nCr(n-1, k-1));
      }
      //u is odd
      else
      {
         k = (u + 1) / 2;
         fu = ((nCr(m-1, k-1) * nCr(n-1, k-2) + nCr(m-1, k-2) * nCr(n-1, k-1)));
      }
      sum += fu;
   }/* end for() */

   p = sum / nCr(m+n, n);

   return p;
}/* end RunsTestCDF() */

/******************************************************************************
nCr()

Computes the number of combinations of n things taken r at a time.
   nCr = n! / [ k! (n-k)! ]
******************************************************************************/
double nCr(int n, int r)
{
   int i;
   double numer, denom;

   numer = denom = 1.00e00;

   if((n - r) < r)
   {
      r = (n - r);
   }

   for(i = 0; i < r; i++) numer *= ((double)(n - i));
   for(i = 0; i < r; i++) denom *= ((double)(r - i));

   return (numer / denom);
}/* end nCr() */

/******************************************************************************
STATS_TestAutocorrelation()

Test out the Durban-Watson and Runs autocorrelation tests.
******************************************************************************/
void STATS_TestAutocorrelation(void)
{
   int i, j;

   //double r[16] = {1,-1,1,1,-1,1,1,1,-1,1,1,1,-1,1,-1,1};
   //double r[20] = {1,-1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,1,1};
    double r[30] = /*{2.86548,1.59472,-0.18592,3.27336,3.98284,1.78204,1.66188,
                  -3.35868,-3.799,-3.52944,-6.62036,-7.18124,-5.352,-2.98244,
                  1.08712,1.18672,2.56656,3.0962,3.45588,6.35548};*/

                  /* Example 2, p-val = 0.025, D = 1.08 
                  {-32.325,-26.598,2.22,-16.962,-1.144,-2.508,-1.963,11.673,
                   -0.509,27.036,-4.419,40.035,23.58,33.943,-2.785,-8.604,
                   0.577,6.849,-18.97,-29.062}; */

                  /* Example 2, p-val = 0.01, D = 0.95 
                  {-27.85108549,-21.94514121,6.992154982,-12.07054883,
                   3.866747361,2.74133974,3.584580214,17.45917259,
                   5.396468783,33.23970926,2.082949731,46.8948383,
                   30.73807877,41.51861544,5.267800196,-0.014366951,
                   9.703465902,16.45265066,-8.829516486,-18.20573935}; */

                   /* example 4, p-val = 0.05, D = 1.08 
                   {1.2911,0.0408,-2.1098,-1.1591,-0.0484,
                    0.2477,0.4824,0.5296,0.9292,0.5075,
                    0.5628,0.4638,-0.2064,-1.5155,-0.0941,
                    0.0000,0.0000,0.0000,0.0000,0.0000}; */

                   /* example 4, p-val = 0.025, D = 0.95 
                  {1.5852,0.3367,-1.7999,-0.8425,0.2645,
                   0.5602,0.8019,0.8535,1.2550,0.8425,
                   0.8908,0.7759,0.1282,-1.1640,0.2173,
                   0.0000,0.0000,0.0000,0.0000,0.0000}; */

                  /* example 4, p-val = 0.01, D = 0.81 
                  {1.7522,0.5047,-1.6238,-0.6626,0.4422,
                  0.7377,0.9834,1.0376,1.4401,1.0329,
                  1.0772,0.9532,0.3184,-0.9641,0.3942,
                  0.0000,0.0000,0.0000,0.0000,0.0000}; */

                  /* SHAZAM example 2, D = 0.57, lwr = 0.00, upr = 1.00 */
                  {-0.166860,-0.141225,-0.056331,-0.097397,-0.218790,
                   -0.180189,0.013461,-0.140583,-0.052185,0.223997,
                   0.246285,0.266461,0.242079,0.073455,0.022708,
                   0.059364,0.097253,0.155865,0.065346,0.073757,
                   0.094268,-0.045548,-0.081874,-0.155034,0.008652,
                  -0.093809,-0.129256,-0.139787,-0.057169,0.113163}; 

                  /* SHAZAM example 1, D = 2.01, lwr = 0.301270, upr = 0.698730 
                  {-5.507930127,-2.576793164,-1.421168928,5.181466085,0.251704188,
                    5.310204304, 1.945793402,-0.574298765,-4.395641673,-1.542614199,
                   -4.59467751,  4.956713133, 8.897180132,-6.415885462, 2.56132361,
                    7.288520277,-9.365001685, 0.00,       0.00,         0.00,
                    0.00,       0.00,         0.00,       0.00,         0.00,
                     0.00,       0.00,         0.00,       0.00,         0.00}; */
                   
   //jacobian/sensitivity matrix
   double J[30][6] = /* Example 1
                      {{1.00,159.3},{1.00,161.2},{1.00,162.8},{1.00,164.6},
                      {1.00,165.9},{1.00,167.9},{1.00,168.3},{1.00,169.7},
                      {1.00,170.5},{1.00,171.6},{1.00,173.9},{1.00,176.1},
                      {1.00,178.0},{1.00,179.1},{1.00,180.2},{1.00,181.2},
                      {1.00,181.6},{1.00,182.5},{1.00,183.3},{1.00,184.3}}; */
                      /* Example 2 
                      {{1.00,75},{1.00,78},{1.00,80},{1.00,82},{1.00,84},
                      {1.00,88},{1.00,93},{1.00,97},{1.00,99},{1.00,104},
                      {1.00,109},{1.00,115},{1.00,120},{1.00,127},{1.00,135},
                     {1.00,144},{1.00,153},{1.00,161},{1.00,170},{1.00,182}}; */
                      /* example 4 
                     {{1.00,794},{1.00,799},{1.00,837},{1.00,855},{1.00,845},
                      {1.00,844},{1.00,863},{1.00,875},{1.00,880},{1.00,905},
                      {1.00,886},{1.00,843},{1.00,904},{1.00,950},{1.00,841},
                      {1.00,0.00},{1.00,0.00},{1.00,0.00},{1.00,0.00},{1.00,0.00}}; */
                     /* SHAZAM example 2 */
                     {{1.00, 25.4, 9.90, 17, 1, 0},
                     {1.00, 26.70, 4.70, 18, 1, 0},
                     {1.00, 29.10, 1.90, 23, 1, 0},
                     {1.00, 29.20, 3.20, 28, 1, 0},
                     {1.00, 29.20, 1.90, 30, 1, 0},
                     {1.00, 27.80, 3.9, 27, 0, 0},
                     {1.00, 27.40, 3.9, 24, 0, 1},
                     {1.00, 28.00, 3.80, 23, 0, 1},
                     {1.00, 28.30, 5.9, 25, 0, 1},
                     {1.00, 28.80, 5.3, 23, 1, 0},
                     {1.00, 29.30, 3.3, 24, 1, 0},
                     {1.00, 29.40, 3.00, 25, 1, 0},
                     {1.00, 29.20, 2.90, 25, 1, 0},
                     {1.00, 29.40, 5.50, 25, 0, 0},
                     {1.00, 30.20, 4.40, 25, 0, 1},
                     {1.00, 31.00, 4.10, 24, 0, 1},
                     {1.00, 31.20, 4.30, 25, 0, 1},
                     {1.00, 31.5, 6.80, 25, 0, 0},
                     {1.00, 31.70, 5.50, 26, 0, 0},
                     {1.00, 32.30, 5.50, 27, 0, 0},
                     {1.00, 32.60, 6.70, 26, 0, 0},
                     {1.00, 32.70, 5.5, 26, 0, 0},
                     {1.00, 33.20, 5.70, 26, 0, 0},
                     {1.00, 33.60, 5.20, 26, 0, 0},
                     {1.00, 34.00, 4.5, 27, 1, 0},
                     {1.00, 34.60, 3.80, 27, 1, 0},
                     {1.00, 35.10, 3.80, 27, 1, 0},
                     {1.00, 35.50, 3.60, 28, 1, 0},
                     {1.00, 36.30, 3.50, 30, 1, 0},
                     {1.00, 36.70, 4.90, 33, 1, 0}}; 
                     /* SHAZAM example 1 
                     {{1,96.7,101},
                     {1,98.1,100.1},
                     {1,100,100},
                     {1,104.9,90.6},
                     {1,104.9,86.5},
                     {1,109.5,89.7},
                     {1,110.8,90.6},
                     {1,112.3,82.8},
                     {1,109.3,70.1},
                     {1,105.3,65.4},
                     {1,101.7,61.3},
                     {1,95.4,62.5},
                     {1,96.4,63.6},
                     {1,97.6,52.6},
                     {1,102.4,59.7},
                     {1,101.6,59.5},
                     {1,103.8,61.3}}; */
                   
   double ** JAC;
   JAC = new double *[30];
   for(i = 0; i < 30; i++) JAC[i] = new double[6];
   for(i = 0; i < 30; i++) for(j = 0; j < 6; j++) JAC[i][j] = J[i][j];

   int pos, neg, runs, clwr, cupr;

   RunsTest(r, 30, &pos, &neg, &runs, &clwr, &cupr);

   for(i = 0; i < 6; i++) delete [] JAC[i];
   delete [] JAC;

}/* end STATS_TestAutocorrelation() */

/******************************************************************************
STATS_TestFdist()

Computes the F-distribution inverse CDF for numerous values of u,v and p.
The results can be compared to those generated in another program (e.g. Excel)
to verify proper operation of the F-distribution code.
******************************************************************************/
void STATS_TestFdist(void)
{
   int i, j, k;
   double val;
   int u[5] = {1,2,4,8,16};
   int v[7] = {1,2,4,8,16,32,64};
   double p[5] ={ 0.250, 0.500, 0.750, 0.975, 0.99};

   printf("***** F-distribution Inverse CDF *****\n");
   printf("u  v  0.25      0.50      0.75      0.975     0.99    \n");
   for(i = 0; i < 5; i++)
   {
      for(j = 0; j < 7; j++)
      {
         printf("%d  %d  ",u[i],v[j]);
         for(k = 0; k < 5; k++)
         {
            val = FdistInvCDF(u[i], v[j], p[k]);
            printf("%lf  ", val);
         }/* end for() */
         printf("\n");
      }/* end for() */
   }/* end for() */
   
}/* end STATS_TestFdist() */

/******************************************************************************
STATS_TestStudentDist()

Computes the t-distribution inverse CDF for numerous values of n,and p.
The results can be compared to those generated in another program (e.g. Excel)
to verify proper operation of the t-distribution code.
******************************************************************************/
void STATS_TestStudentDist(void)
{
   int i, j;
   double val;
   int u[11] = {1,2,4,8,16,32,64,128,256,512,1024};
   double p[5] ={ 0.550, 0.750, 0.800, 0.975, 0.995};

   printf("***** Student's t-distribution Inverse CDF *****\n");
   printf("u  0.550      0.750      0.800      0.975     0.995\n");
   for(i = 0; i < 11; i++)
   {
      printf("%d  ",u[i]);
      for(j = 0; j < 5; j++)
      {         
         val = StudentInvCDF(u[i], p[j]);
         printf("%lf  ", val);
      }/* end for() */
      printf("\n");
   }/* end for() */
}/* end STATS_TestStudentDist() */

/******************************************************************************
STATS_TestStdNormDist()

Computes the std. normal distribution inverse CDF for numerous values of p.
The results can be compared to those generated in another program (e.g. Excel)
to verify proper operation of the code.
******************************************************************************/
void STATS_TestStdNormDist(void)
{
   int i;
   double val;   

   printf("***** Standard Normal Inverse CDF *****\n");
   printf("Probability  X\n");
   for(i = 0; i < 19; i++)
   {
      printf("%lf  ",0.05*(i+1));
      val = StdNormInvCDF(0.05*(i+1));
      printf("%lf\n", val);
   }/* end for() */

   printf("%lf  ",0.99999);
   val = StdNormInvCDF(0.99999);
   printf("%lf\n", val);
}/* end STATS_TestStdNormDist() */
