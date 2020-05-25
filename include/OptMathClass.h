/******************************************************************************
File     : OptMathClass.h
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The OptMathClass is used to compute mathematical measures, namely 1st and 2nd 
order derivatives, that are used in certain optimization algorithms.

Version History
08-22-03    lsm   created file
03-24-04    lsm   Added CalcOptimalStepSize() option
08-18-04    lsm   Added metrics collection and reporting, memory fragmentation
                  fixes
03-21-05    lsm   Added support for parameter-specific relative increments
01-01-07    lsm   Added support for additional FD increment types. OptMathClass
                  now uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef OPT_MATH_CLASS_H
#define OPT_MATH_CLASS_H

#include "MyHeaderInc.h"

// forward decs
class ModelABC;

/******************************************************************************
class OptMathClass

******************************************************************************/
class OptMathClass
{      
   public:
      OptMathClass(ModelABC * pModel);
      ~OptMathClass(void){ DBG_PRINT("OptMathClass::DTOR"); Destroy(); }
      void Destroy(void);

      Unchangeable1DArray CalcGradient(double * fmin, double * pmin);
      Unchangeable2DArray CalcHessian(void);
      void WriteMetrics(FILE * pFile);

   private:
      //Configuration variables (to be input by user)
      FiniteDiffType m_DiffType; //finite diff. method used to calc. derivatives     
      FiniteDiffIncType m_DiffIncType; //increment type used to calc. derivatives     
      double * m_pDiffInc; //finite diff. increment array
      double   m_MinInc; //smallest allowable increment
      
      double *  m_pGrad; // gradient
      double ** m_pHess; // hessian                       

      //points in the design space, used by various routines
      double * m_pHessPoint; 
      double * m_pGradPoint; 
      double * m_pStepPoint; 
      double * m_pDiffPoint;
            
      int m_NumParams;  //number of parameters

      ModelABC * m_pModel;

      //metrics
      int m_DiffCount;
      int m_GradCount;
      int m_StepCount;
      int m_HessCount;

      void InitFromFile(IroncladString pMathFileName); 

      /* calculates derivative at the current model point with 
      respect to the parameter located at parmIdx */
      double CalcDerivative(int parmIdx, double * fmin, double * pmin);

      /* calculates optimal step size at the current model 
      point with respect to the parameter located at parmIdx */
      double CalcOptimalStepSize(int idx);
}; /* end class OptMathClass */

extern "C" {
void Gradient_Program(int argc, StringType argv[]);
void Hessian_Program(int argc, StringType argv[]);
}

#endif /* OPT_MATH_CLASS_H */

