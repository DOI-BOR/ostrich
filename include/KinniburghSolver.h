/******************************************************************************
File     : KinniburghSolver.h
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The KinniburghSolver class solves the non-linear Kinniburgh isotherm equation.
See Kinniburgh.cpp for details.

Version History
07-28-07    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef KINNIBURGH_SOLVER_H
#define KINNIBURGH_SOLVER_H

#include "MyHeaderInc.h"

//forward declarations
class IsothermABC;
class ObservationGroup;

/******************************************************************************
class KinniburghSolver
******************************************************************************/
class KinniburghSolver
{   
   public:
      KinniburghSolver(IsothermABC * pIso);
     ~KinniburghSolver(void){ DBG_PRINT("KinniburghSolver::DTOR"); Destroy(); }
      void Destroy(void);
      void Compute(void);
      void Compute(ObservationGroup * pObs);
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
      double m_Cupr;
      double m_Clwr;
      IsothermABC * m_pIso;
}; /* end class KinniburghSolver */

#endif /* KINNIBURGH_SOLVER_H */


