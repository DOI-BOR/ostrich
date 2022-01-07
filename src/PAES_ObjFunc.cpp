/******************************************************************************
File       : PAES_ObjFunc.cpp
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)

PAES - Pareto archived evolution strategy

Version History
12-29-17    lsm   created file
******************************************************************************/

#include "PAES_ObjFunc.h"
#include "PAES_Genetics.h"
#include "string.h"
#include "ParameterGroup.h"
#include "ObjectiveFunction.h"
#include "ParameterABC.h"

/******************************************************************************
PAES_MultiobjectiveProblem

CTOR
******************************************************************************/
PAES_MultiobjectiveProblem::PAES_MultiobjectiveProblem(void)
{
   m_variable = NULL;
} /* end CTOR */

/******************************************************************************
PAES_MultiobjectiveProblem::adjustPrecision
******************************************************************************/
void PAES_MultiobjectiveProblem::adjustPrecision(int variable)
{
   double U,L, P;
   int representationBase = 2;
   int bitsPerVar;
   U = m_upperLimit[variable];
   L = m_lowerLimit[variable];
   P = m_precision[variable];
   if(U <= L)
   {
      bitsPerVar = 1;
   }
   else
   {
      bitsPerVar = (int)ceil(log((U - L) * pow(10.0, P) * 1.000001) / log((double)representationBase));
   }
   if(bitsPerVar < 1)
   {
      bitsPerVar = 1;
   }
   m_bitsPerVariable[variable] = bitsPerVar;
} /* end MultiobjectiveProblem::adjustPrecision() */

/******************************************************************************
PAES_MultiobjectiveProblem::constraintsAreSatisfied
******************************************************************************/
bool PAES_MultiobjectiveProblem::constraintsAreSatisfied(PAES_Individual * individual)
{
   return true;
} /* end MultiobjectiveProblem::constraintsAreSatisfied() */

/******************************************************************************
PAES_MultiobjectiveProblem::numberOfNonSatisfiedConstraints
******************************************************************************/
int PAES_MultiobjectiveProblem::numberOfNonSatisfiedConstraints(PAES_Individual * ind)
{
   return 0; 
} /* end  MultiobjectiveProblem::numberOfNonSatisfiedContraints() */

/******************************************************************************
PAES_MultiobjectiveProblem::initializeRealVariableType
******************************************************************************/
void PAES_MultiobjectiveProblem::initializeRealVariableType(PAES_VariableType variableType)
{
   int i;
   for (i = 0; i < m_numberOfVariables; i++)
   {
      if (variableType == PAES_REAL) 
      {
         m_variable[i] = PAES_REAL;
      }
      else if (variableType == PAES_BINARY_REAL) 
      {
         m_variable[i] = PAES_BINARY_REAL;
         adjustPrecision(i);
         fprintf(m_cout, " bits variable %d : %d\n", i, m_bitsPerVariable[i]);
      } /* end else if() */
      else 
      {
         fprintf(m_cout, "MultiobjectiveProblem::initializeVariableType->\n");
         fprintf(m_cout, "   The type of the variable must be Real or BinaryReal.\n");
         fprintf(m_cout, "   Type received: %d\n", variableType);
         fclose(m_cout);
         exit(-1);
      } /* end else() */
   }/* end for() */
} /* end MultiobjectiveProblem::initializeVariableType() */

/******************************************************************************
PAES_MultiobjectiveProblem::print

Prints the information associated to a multiobjective problem
******************************************************************************/
void PAES_MultiobjectiveProblem::print(void) 
{
   int i;

   fprintf(m_cout, "Problem Name      : %s\n", m_problemName);
   fprintf(m_cout, "Variables         : %d\n", m_numberOfVariables);
   fprintf(m_cout, "Functions         : %d\n", m_numberOfFunctions);
   fprintf(m_cout, "\n");

   fprintf(m_cout, "Limits            : ");
   for (i = 0; i < m_numberOfVariables; i++) 
   {
      fprintf(m_cout, "(%E , %E)\n", m_lowerLimit[i], m_upperLimit[i]);
   } /* end for() */
   fprintf(m_cout, "\n");

   fprintf(m_cout, "====\n");
} /* end MultiobjectiveProblem::print() */

/******************************************************************************
Destroy()
******************************************************************************/
void PAES_OstrichIf::Destroy(void)
{
   delete [] m_variable;
   delete [] m_upperLimit;
   delete [] m_lowerLimit;
   delete [] m_precision;
   delete [] m_bitsPerVariable;
}/* end PAES_OstrichIf::Destroy() */

/******************************************************************************
PAES_OstrichIf

CTOR
******************************************************************************/
PAES_OstrichIf::PAES_OstrichIf(PAES_VariableType variableType, ModelABC * pModel, FILE * pOut)
{
   ParameterGroup * pGroup;
   int i;

   m_cout = pOut;
   strcpy(m_problemName, "OstrichIf");

   m_pModel = pModel;
   pGroup = pModel->GetParamGroupPtr();

   m_numberOfVariables = pGroup->GetNumParams();
   m_numberOfFunctions = pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);
   m_numberOfConstraints = 0;
    
   m_upperLimit = new double[m_numberOfVariables];
   m_lowerLimit = new double[m_numberOfVariables];
   m_precision = new int[m_numberOfVariables];
   m_bitsPerVariable = new int[m_numberOfVariables];

   for(i = 0; i < m_numberOfVariables; i++)
   {
      m_lowerLimit[i] = pGroup->GetParamPtr(i)->GetLwrBnd();
      m_upperLimit[i] = pGroup->GetParamPtr(i)->GetUprBnd();
      m_precision[i] = 5;
      m_bitsPerVariable[i] = 8;
   }/* end for() */
   
   m_variable = new PAES_VariableType[m_numberOfVariables];

   initializeRealVariableType(variableType);
   fprintf(m_cout, "Created a %s problem\n", m_problemName);
} /* end CTOR */

/******************************************************************************
PAES_OstrichIf::evaluate()
******************************************************************************/
void PAES_OstrichIf::evaluate(PAES_Individual *individual, double * X, int np, double * F, int nf) 
{
   int i;

   /* retrieve variables from gene encoding */
   for (i = 0; i < m_numberOfVariables; i++) 
   {
      X[i] = (individual->m_chromosome->m_gene[i])->getRealAllele();
   }/* end for() */

   /* pass variables along to model */
   m_pModel->GetParamGroupPtr()->WriteParams(X);

   /* run model */
   m_pModel->Execute(F, nf);

   /* pass along results */
   for (i = 0; i < nf; i++)
   {
      individual->m_fitness[i] = F[i];
   }/* end for() */
} /* PAES_OstrichIf::evaluate() */


