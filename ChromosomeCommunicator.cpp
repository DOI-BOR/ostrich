/******************************************************************************
File     : ChromosomeCommunicator.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

Because Chromosomes are a coded version of design variables, translation is 
necessary when information is exchanged between a Genetic Algorithm and a 
Model. Therefore, the ModelChromoComm class acts as an interface 
between the Model and Chromosome classes.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field, updated comments and renamed
                  the class.
07-08-04    lsm   switched over to using ParameterABC class
08-17-04    lsm   RAM fragmentation fixes
10-19-04    lsm   Added support for binary coded GA
01-01-07    lsm   Algorithm uses abstract model base class (ModelABC).
******************************************************************************/
#include <stdio.h>
#include <math.h>

#include "ChromosomeCommunicator.h"
#include "Chromosome.h"
#include "ModelABC.h"
#include "ModelBackup.h"
#include "Gene.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
GetParamGroupPtr()
******************************************************************************/
ParameterGroup * ModelChromoComm::GetParamGroupPtr(void)
{
	return  m_pModel->GetParamGroupPtr(); 
}/* end GetParamGroupPtr() */

/******************************************************************************
CTOR

Assigns model pointer.
******************************************************************************/
ModelChromoComm::ModelChromoComm(ModelABC * pModel)
{
   m_pModel = pModel;
   m_xb = NULL;
   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()
******************************************************************************/
void ModelChromoComm::Destroy(void)
{
   delete [] m_xb;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
EvalFitness()

Evaluates the fitness of a chromosme and stores the result in the chromosomes
fitness variable.
******************************************************************************/
void ModelChromoComm::EvalFitness(Chromosome * pChromo)
{   
   static double fbest = HUGE_VAL;
   double fitness, val;
   Gene * pGene;
   ParameterGroup * pParamGroup;
   ParameterABC * pParam;
   int i, numGenes, numParams;

   pParamGroup = m_pModel->GetParamGroupPtr();
   numGenes = pChromo->GetNumGenes();
   numParams = pParamGroup->GetNumParams();

  if(m_xb == NULL)
  {
     m_xb = new double[numParams];
     pParamGroup->ReadParams(m_xb);
  }

   //sanity check (should have same number of genes as parameters)
   if(numParams != numGenes)
   {
      LogError(ERR_MISMATCH, "Number of genes != Number of parameters");
      ExitProgram(1);
   }/* end if() */

   for(i = 0; i < numGenes; i++)
   {      
      pGene = pChromo->GetGenePtr(i);
      pParam = pParamGroup->GetParamPtr(i);
      val = pGene->GetValue();
      pParam->SetEstVal(val);      
   }/* end for() */

   fitness = -1.00 * m_pModel->Execute();

   if((-fitness) <= fbest)
   {
      fbest = -fitness;
      pParamGroup->ReadParams(m_xb);
   }
   pChromo->SetFitness(fitness);
} /* end EvalFitness() */

/******************************************************************************
MakeParameterCorrections()

Corrects the genes in a chromosme according to expert rules.
******************************************************************************/
void ModelChromoComm::MakeParameterCorrections(Chromosome * pChromo)
{   
   static double a = 0.00;
   double val, lwr, upr;
   Gene * pGene;
   ParameterGroup * pParamGroup;
   ParameterABC * pParam;
   int i, numGenes, numParams;

   pParamGroup = m_pModel->GetParamGroupPtr();
   numGenes = pChromo->GetNumGenes();
   numParams = pParamGroup->GetNumParams();

  if(m_xb == NULL)
  {
     m_xb = new double[numParams];
     pParamGroup->ReadParams(m_xb);
  }

   //sanity check (should have same number of genes as parameters)
   if(numParams != numGenes)
   {
      LogError(ERR_MISMATCH, "Number of genes != Number of parameters");
      ExitProgram(1);
   }/* end if() */

   for(i = 0; i < numGenes; i++)
   {      
      pGene = pChromo->GetGenePtr(i);
      pParam = pParamGroup->GetParamPtr(i);
      lwr=pParam->GetLwrBnd();
      upr=pParam->GetUprBnd();
      val = pGene->GetValue();
      val = TelescopicCorrection(lwr, upr, m_xb[i], a, val);
      pGene->SetValue(val);
      pParam->SetEstVal(val);      
   }/* end for() */

   //inerface with expert judgement module
   m_pModel->PerformParameterCorrections();
   for(i = 0; i < numGenes; i++)
   {
      pGene = pChromo->GetGenePtr(i);
      pParam = pParamGroup->GetParamPtr(i);
      val = pParam->GetEstVal();
      pGene->SetValue(val);
   }/* end for() */

   a += 1.00/(double)m_MaxEvals;
} /* end MakeParameterCorrections() */

/******************************************************************************
CreateProto()

Creates a prototype chromosome using the information from the models paramter 
group.
******************************************************************************/
Chromosome * ModelChromoComm::CreateProto(double rate)
{ 
   Chromosome * pChromo;
   Gene * pGene;
   ParameterGroup * pParamGroup;   
   ParameterABC * pParam;
   int i, numParams;
   double upr,lwr, xover;

   pParamGroup = m_pModel->GetParamGroupPtr();
   numParams = pParamGroup->GetNumParams();

   // set gene xover rate to achieve a 
   // net crossover of 0.5 per chomosome
   xover=1-(pow(0.5,(1.0/(double)numParams)));

   NEW_PRINT("Chromosome", 1);
   pChromo = new Chromosome(0.00, numParams);
   MEM_CHECK(pChromo);
   
   for(i = 0; i < numParams; i++)
   {                  
      pParam = pParamGroup->GetParamPtr(i);
      lwr = pParam->GetLwrBnd();
      upr = pParam->GetUprBnd();

      if(GetProgramType() == GA_PROGRAM) //real encoded genes
      {
         NEW_PRINT("RealEncodedGene", 1);
         pGene = new RealEncodedGene(0.5*(upr+lwr), lwr, upr, rate, xover);
         MEM_CHECK(pGene);
      }
      else //binary encoded genes
      {
         NEW_PRINT("BinaryEncodedGene", 1);
         pGene = new BinaryEncodedGene(0.5*(upr+lwr), lwr, upr, rate, xover);
         MEM_CHECK(pGene);
      }

      pChromo->SetGenePtr(pGene, i);
   }/* end for() */

   return pChromo;
}/* end CreateProto() */

/******************************************************************************
ConvertChromosome()

Uses the supplied chromosome to return the resulting parameter group.
******************************************************************************/
ParameterGroup * ModelChromoComm::ConvertChromosome(Chromosome * pChromo)
{
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   int numGenes;
   int numParams;
   int i;
   double val;
   Gene * pGene;

   pGroup = m_pModel->GetParamGroupPtr();
   numGenes = pChromo->GetNumGenes();
   numParams = pGroup->GetNumParams();

   if(numParams != numGenes)
   {
      LogError(ERR_MISMATCH, "Number of genes != Number of parameters");
      ExitProgram(1);
   }/* end if() */

   for(i = 0; i < numGenes; i++)
   {
      pGene = pChromo->GetGenePtr(i);      
      val = pGene->GetValue();
      pParam = pGroup->GetParamPtr(i);
      pParam->SetEstVal(val);
   }/* end for() */

   return pGroup;
}/* end ConvertChromosome() */


