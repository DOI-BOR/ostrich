/******************************************************************************
File     : McCammonSolver.h
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The McCammonSolver class solves the non-linear McCammon equation for handling
errors in both variables (q and C).

See McCammon.cpp for details.

Version History
07-28-07    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef MCCAMMON_SOLVER_H
#define MCCAMMON_SOLVER_H

#include "MyHeaderInc.h"

//forward decs
class IsothermABC;
class ObservationGroup;

/******************************************************************************
class McCammonSolver
******************************************************************************/
class McCammonSolver
{   
   public:
      McCammonSolver(IsothermABC * pIso);
     ~McCammonSolver(void){ DBG_PRINT("McCammonSolver::DTOR"); Destroy();}
      void Destroy(void);
      void Compute(void);
      void Compute(ObservationGroup * pObs);
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
}; /* end class McCammonSolver */

#endif /* MCCAMMON_SOLVER_H */


