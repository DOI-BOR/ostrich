/******************************************************************************
File       : PAES_ObjFunc.h
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)

PAES - Pareto archived evolution strategy

Version History
12-29-17    lsm   created file
******************************************************************************/

#ifndef PAES_OBJ_FUNC_H
#define PAES_OBJ_FUNC_H

#include "PAES_Helpers.h"
#include "ModelABC.h"

//forward declarations
class PAES_Individual;

/******************************************************************************
class PAES_MultiobjectiveProblem

A base class for objective functions to be solved by PAES
******************************************************************************/
class PAES_MultiobjectiveProblem 
{
   public:
      int m_numberOfVariables;
      int m_numberOfFunctions;
      int m_numberOfConstraints;

      char m_problemName[DEF_STR_SZ];
  
      PAES_VariableType * m_variable; //!< Types of the decision variables
  
      double * m_upperLimit; //!< To be used with continuous variables
      double * m_lowerLimit; //!< To be used with continuous variables
  
      int    * m_precision; //!< To be used by binary variables
      int      m_numberOfBits; //!< To be used by binary variables
      int    * m_bitsPerVariable; //!< To be used by binary variables

      FILE * m_cout;
  
      PAES_MultiobjectiveProblem(void);
    
      virtual void evaluate(PAES_Individual * individual, double * X, int np, double * F, int nf) = 0;
      virtual bool constraintsAreSatisfied(PAES_Individual * individual);
      virtual int  numberOfNonSatisfiedConstraints(PAES_Individual * individual);
  
      // Methods for binary variables
      void adjustPrecision(int variable);

      // Methods for real variables
      void initializeRealVariableType(PAES_VariableType variableType);
  
      void print(void);
}; /* end class PAES_MultiobjectiveProblem */

/******************************************************************************
class PAES_OstrichIf
 
The interface between OSTRICH objects and PAES objects
 ******************************************************************************/
class PAES_OstrichIf : public PAES_MultiobjectiveProblem
{
   public:
      // Constructor
      PAES_OstrichIf(PAES_VariableType type, ModelABC * pModel, FILE * pOut);
      ~PAES_OstrichIf(void){ DBG_PRINT("PAES_OstrichIf::DTOR"); Destroy(); }
  
      void evaluate(PAES_Individual *individual, double * X, int np, double * F, int nf);
      void Destroy(void);

   private:
      ModelABC * m_pModel;
}; /* end class PAES_OstrichIf */

#endif /* PAES_OBJ_FUNC_H */
