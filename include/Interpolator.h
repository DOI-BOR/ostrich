/******************************************************************************
File      : Interpolator.h
Author    : L. Shawn Matott and James Craig
Copyright : 2004, L. Shawn Matott and James Craig

A class that uses radial basis function to interpolates the value at a given 
point, based on known values at a finite set of points.

Version History
05-25-04    lsm   created
******************************************************************************/
#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "MyHeaderInc.h"

/******************************************************************************
class Interpolator
******************************************************************************/
class Interpolator
{
   public:
      Interpolator(int nmax);
      ~Interpolator(void){ DBG_PRINT("Interpolator::DTOR"); Destroy(); }
      void SetBasis(MyPoint * vals, int n);
      double Evaluate(MyPoint * point);
      double Interpolate(bool debug);
      void Destroy(void);

   private:      
      MyPoint * m_pBasis; //known values
      double ** m_A; //lhs matrix (A)
      double ** m_Ainv; //inverse of A
      double *  m_b; //rhs vector (b)
      int m_Order; //number of known values and coeffs.
      int m_MaxOrder;//maximum number of known values used by the interpolator
      double * m_pCoeffs; //coefficients
      double * m_pRadius; //radii of MQ-RBF expansion
      double m_AvgVal; //average z value form list of known values
}; /* end class Interpolator */

#endif /* INTERPOLATOR_H */



