/******************************************************************************
File      : TelescopingBounds.cpp
Author    : L. Shawn Matott
Copyright : 2012, L. Shawn Matott

Encapsulates a group of functions for dynamically adjusting parameter bounds.

Version History
12-12-12    lsm   added copyright information and initial comments.
******************************************************************************/
#include <math.h>

#include "TelescopingBounds.h"

#include "MyTypes.h"

double flin(double a, double b)
{
  return 1.0-0.99*a;
}

double fpvx(double a, double b)
{
  return pow(10.00, -2.00*a);
}

double fdcv(double a, double b)
{
  b = 0.2;
  if (a <= b) return 1.00;
  return 0.01+0.99*sin(MY_PI*0.5*(1.00 - a)/(1.00 - b));
}

double fcve(double a, double b)
{
  return 0.01+0.99*sin(MY_PI*0.5*(1.00 - a));
}

double fvex(double a, double b)
{
  return 1.00+0.99*sin(MY_PI*0.5*(4.00 - a));
}

double telescope_parameter(double xmin, double xmax, double xbest, 
                         double a, double xnew, double (*f)(double, double))
{
  double xmin_new, xmax_new, xp, range;
  revise_bounds(xmin, xmax, xbest, a, &xmin_new, &xmax_new, f);
  if(( xnew > xmin_new) && (xnew < xmax_new))
    return xnew;

  xp = (xnew - xmin)/(xmax - xmin);
  range = (xmax_new - xmin_new);
  xnew = xmin_new + xp*range;
  return xnew;
}

void revise_bounds(double xmin, double xmax, double xbest, 
                   double a, double * xmin_new, double * xmax_new, double (*f)(double, double))
{
  double diff;
  double range = xmax - xmin;
  if(a < 0.00) a = 0.00;
  if(a > 1.00) a = 1.00;
  *xmin_new = xbest - f(a,0)*range*(0.5);
  *xmax_new = xbest + f(a,0)*range*(0.5);

   if(*xmin_new < xmin) 
   {
     diff = xmin - *xmin_new;
     *xmax_new += diff;
     *xmin_new = xmin;
   }

    if(*xmax_new > xmax) 
    {
      diff = *xmax_new - xmax;
      *xmin_new -= diff;
      *xmax_new = xmax;
    }    
}


