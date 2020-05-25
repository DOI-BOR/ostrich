/******************************************************************************
File     : OptSearchClass.h
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The OptSearchClass is used to perform one-dimensional searches for the optimum
step size of a given search direction.

Version History
08-25-03    lsm   created file
08-18-04    lsm   Added metrics collection and reporting, memory fragmentation
                  fixes, Golden Section search modifications
01-01-07    lsm   OptSearchClass now uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef OPT_SEARCH_CLASS_H
#define OPT_SEARCH_CLASS_H

#include "MyHeaderInc.h"

// forward decs
class ModelABC;

typedef struct MIN_BRACKET_STRUCT
{
   double a, b, c;
   double fa, fb, fc;
}MinBracketStruct;

typedef struct MIN_PT_STRUCT
{
   double x, fx;
}MinPtStruct;

/* User can choose from these 1-D search methods */
#define GSECT_SEARCH (0)
#define BRENT_SEARCH (1)

/******************************************************************************
class OptSearchClass

******************************************************************************/
class OptSearchClass
{      
   public:     
      OptSearchClass(ModelABC * pModel);
      ~OptSearchClass(void){ DBG_PRINT("OptSearchClass::DTOR"); Destroy(); }
      void Destroy(void);
      double CalcStepSize(Unchangeable1DArray pDir, double * fmin, double * xmin);
      void WriteMetrics(FILE * pFile);

   private:
      MinPtStruct * GoldSect(MinBracketStruct * pBrack, double * minf, double * minp);
      MinPtStruct * Brent(MinBracketStruct * pBrack, double * minf, double * minp);
      MinBracketStruct * BracketMinimum(double a, double b, double * fmin, double * xmin);
      double CalcF(double alpha,double * fmin, double * xmin);
      void InitFromFile(IroncladString pFileName);
      double LimitStepSize(double alpha, double * pDir);

      ModelABC * m_pModel;

      //frequently used data structs
      double           * m_pStepPoint; //CalcStepSize() point
      double           * m_pAlphaPoint; //CalcF() points
      double           * m_pStartPoint; 
      MinBracketStruct * m_pMinBrack;
      MinPtStruct      * m_pMinPt;

      int m_NumParams;  // number of parameters
      double * m_pDir;  // search direction
      double m_Step;    // step size     

      // convergence criteria for 1-D search
      double m_SearchConvVal;
      //type of 1-D search
      int m_SearchType;

      //metrics
      int m_BoundMinCount;
      int m_GoldSectCount;
      int m_BrentCount;
}; /* end class OptSearchClass */

#endif /* OPT_SEARCH_CLASS_H */

