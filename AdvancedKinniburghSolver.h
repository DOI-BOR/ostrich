/******************************************************************************
File     : AdvancedKinniburghSolver.h
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The AdvancedKinniburghSolver class solves the non-linear Kinniburgh isotherm equation.
See AdvancedKinniburgh.cpp for details.

Version History
03-10-10    lsm   created from KinniburghSolver.h
******************************************************************************/
#ifndef ADV_KINNIBURGH_SOLVER_H
#define ADV_KINNIBURGH_SOLVER_H

#include "MyHeaderInc.h"

//forward class delcarations
class IsothermABC;

/******************************************************************************
class AdvancedKinniburghSolver
******************************************************************************/
class AdvancedKinniburghSolver
{   
   public:
      AdvancedKinniburghSolver(IsothermABC * pIso, double X);
     ~AdvancedKinniburghSolver(void){ DBG_PRINT("AdvancedKinniburghSolver::DTOR"); Destroy();}
      void Destroy(void);
      void Compute(void);
      bool Initialize(char * pStr);
  
   private:
      double BisectionSearch(int i);
      double F(double C, double A, double BD);

      int m_NumOut;
      int m_MaxIters; //maximum number of bisections
      char * m_pOutFile;
      double * m_pC;  //array of aqueous concentrations
      double * m_pA; //experimental constant
      double * m_pB; //experimental constant
      double * m_pD; //experimental constant
      double m_X;    //loss term
      double m_Cupr;
      double m_Clwr;
      IsothermABC * m_pIso;
}; /* end class AdvancedKinniburghSolver */

#endif /* ADV_KINNIBURGH_SOLVER_H */


