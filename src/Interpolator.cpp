/******************************************************************************
File      : Interpolator.h
Author    : L. Shawn Matott and James Craig
Copyright : 2004, L. Shawn Matott and James Craig

A class that uses radial basis function to interpolates the value at a given 
point, based on known values at a finite set of points.

Version History
05-25-04    lsm   created
******************************************************************************/
#include <math.h>
#include <string.h>

#include "Interpolator.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
CTOR
******************************************************************************/
Interpolator::Interpolator(int nmax)
{
   int i;
   m_MaxOrder = nmax;

   //allocate memory
   m_b = new double[nmax];
   m_A = new double*[nmax];
   m_Ainv = new double*[nmax];
   for(i = 0; i < nmax; i++)
   {
      m_A[i] = new double[nmax];
      m_Ainv[i] = new double[nmax];
   }
   MEM_CHECK(m_Ainv[i-1]);
   m_pCoeffs = new double[nmax];
   m_pRadius = new double[nmax];
   MEM_CHECK(m_pRadius);
} /* end CTOR() */

/******************************************************************************
Interpolate()

Compute the MQ-RBF (multiquadric radial basis function) coefficients.
******************************************************************************/
double Interpolator::Interpolate(bool bDebug)
{
   double * vi, * vj, rj, dij, gij;
   int i, j, k, ndim;
   bool test;

   /*------------------------------------------
   Set up system of linear equations (Ax = b) 
   and solve.   
   ------------------------------------------*/

   //assign A and b
   ndim = m_pBasis[0].ndim;
   for(i = 0; i < m_Order; i++)
   {
      vi = m_pBasis[i].v;

      m_b[i] = m_pBasis[i].F - m_AvgVal;
      for(j = 0; j < m_Order; j++)
      {
         vj = m_pBasis[j].v;
         rj = m_pRadius[j];
         dij = 0.00;
         for(k = 0; k < ndim; k++)
         {
            dij += ((vi[k] - vj[k])*(vi[k] - vj[k]));
         }
         gij = sqrt(dij + (rj*rj));

         m_A[i][j] = gij;
         m_Ainv[i][j] = 0.00;         
      }/* end for() */
   }/* end for() */

   //invert A
   if(bDebug == true) fprintf(stdout, "Inverting matrix\n");
   test = MatInv(m_A, m_Ainv, m_Order);

   if(test == false)
   {
      printf("unable to interpolate.\n");
      return 0.00;
   }

   //pre-multiply b by Ainv and assign result to coeff. array
   if(bDebug == true) fprintf(stdout, "Computing coefficients\n");
   VectMult(m_Ainv, m_b, m_pCoeffs, m_Order, m_Order);   
         
   return 0.00;
} /* end Interpolate() */

/******************************************************************************
Evaluate()

Evualute z = f(x1,x2,x3,...,xn) using the MQ-RBF interpolation.
******************************************************************************/
double Interpolator::Evaluate(MyPoint * point)
{
   int i, j, ndim;
   double sum, ai, gi, di, ri;
   double * v, *vi;

   sum = 0.00;
   v = point->v;
   ndim = point->ndim;
   for(i = 0; i < m_Order; i++)
   {
      ai = m_pCoeffs[i];
      ri = m_pRadius[i];
      vi = m_pBasis[i].v;
      di = 0.00;
      for(j = 0; j < ndim; j++)
      {
         di += ((v[j] - vi[j])*(v[j] - vi[j]));
      }
      gi = sqrt(di + (ri*ri));

      sum += (ai * gi);
   }/* end for() */

   sum += m_AvgVal;

   point->F = sum;

   return sum;
} /* end Evaluate() */

/******************************************************************************
Destroy

Destructor.
******************************************************************************/
void Interpolator::Destroy(void)
{
   int i;
   delete [] m_pCoeffs;
   delete [] m_pRadius;

   /*------------------------------------------
   Free up matrices and vectors
   ------------------------------------------*/   
   delete [] m_b;
   for(i = 0; i < m_MaxOrder; i++)
   {
      delete [] m_A[i];
      delete [] m_Ainv[i];
   }
   delete [] m_A;
   delete [] m_Ainv;
} /* end DTOR */

/******************************************************************************
SetBasis()

Assign a new set of known points to interpolate from.
******************************************************************************/
void Interpolator::SetBasis(MyPoint * vals, int n)
{
   int i;
   m_pBasis = vals;
   m_Order = n;
   if(n > m_MaxOrder) n = m_MaxOrder;

   //compute average F value
   m_AvgVal = 0.00;
   for(i = 0; i < n; i++){ m_AvgVal += vals[i].F;}
   m_AvgVal /= (double)n;

   //zero coeffs.
   for(i = 0; i < n; i++){ m_pCoeffs[i] = 0.00;}

   //zero radii
   for(i = 0; i < n; i++){ m_pRadius[i] = 0.00;}
}/* end SetBasis() */
