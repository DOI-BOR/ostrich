/******************************************************************************
File     : OrearSolver.h
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The OrearSolver class solves the non-linear Orear equation for handling
errors in both variables (q and C).

See Orear.cpp for details.

Version History
07-28-07    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef OREAR_SOLVER_H
#define OREAR_SOLVER_H

#include "MyHeaderInc.h"

// forward decs
class IsothermABC;

/******************************************************************************
class OrearSolver
******************************************************************************/
class OrearSolver
{   
   public:
      OrearSolver(IsothermABC * pIso);
     ~OrearSolver(void){ DBG_PRINT("OrearSolver::DTOR"); Destroy(); }
      void Destroy(void);
      void Compute(void);
      bool Initialize(char * pStr);
  
   private:
      double BisectionSearch(int i);
      double F(double C, double Cobs, double qobs, double wc, double wq);

      int m_NumOut;
      int m_MaxIters; //maximum number of bisections
      char * m_pOutFile;
      double * m_pC;  //array of aqueous concentrations
      double * m_pq; //array of sorbed concentrations
      double * m_pWc;  //array of aqueous obs. weights
      double * m_pWq; //array of sorbed obs. weights
      double m_Cupr;
      double m_Clwr;
      IsothermABC * m_pIso;
}; /* end class OrearSolver */

#endif /* OREAR_SOLVER_H */


