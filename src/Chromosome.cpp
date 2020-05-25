/******************************************************************************
File     : Chromosome.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

A Chromosome is a set of design variables (also called genes) that make up a 
single solution to a given optimization problem. As their name implies, 
chromosomes are used by the Genetic Algorithm.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes
******************************************************************************/
#include "Chromosome.h"
#include "Gene.h"

#include "Exception.h"

/******************************************************************************
CTOR

Assigns member variables and allocates the initial Gene List.
******************************************************************************/
Chromosome::Chromosome(double fitness, int numGenes)
{
   int i;
   
   m_NumGenes = numGenes;
   m_Fitness  = fitness;   
   NEW_PRINT("Gene *", numGenes);
   m_pGenes   = new Gene * [numGenes];
   MEM_CHECK(m_pGenes);

   for(i = 0; i < numGenes; i++) { m_pGenes[i] = NULL; }
   IncCtorCount();
}/* end CTOR */

/******************************************************************************
SetFitness()

Assigns a fitness value.
******************************************************************************/
void Chromosome::SetFitness(double fitness)
{
  m_Fitness = fitness;
}/* end SetFitness() */

/******************************************************************************
GetFitness()

Retrieves the fitness value.
******************************************************************************/
double Chromosome::GetFitness(void)
{
  return m_Fitness;
}/* end GetFitness() */

/******************************************************************************
GetGenePtr()

Retrieves the Gene located at ith index of the Gene List. If the index is out 
of bounds, the program will be aborted.
******************************************************************************/
Gene * Chromosome::GetGenePtr(int i)
{
   //sanity check on index
   if(i >= m_NumGenes) 
   {
      LogError(ERR_ARR_BNDS, "GetGenePtr(): index out of bounds");
      ExitProgram(1);
   }/* end if() */

   return m_pGenes[i];
}/* end GetGenePtr() */

/******************************************************************************
SetGenePtr()

Assigns the Gene located at ith index of the Gene List. If the index is out of 
bounds, the program will be aborted. If the index is already assigned to a 
Gene, the old Gene will freed.
******************************************************************************/
void Chromosome::SetGenePtr(Gene * pGene, int i)
{
   //sanity check on index
   if(i >= m_NumGenes) 
   {
      LogError(ERR_ARR_BNDS, "SetGenePtr(): index out of bounds");
   }/* end if() */
   else
   {
      delete m_pGenes[i];
      m_pGenes[i] = pGene;
   }   
}/* end SetGenePtr() */

/******************************************************************************
GetNumGenes()

Retrieves the number of genes in the Gene List. 
******************************************************************************/
int Chromosome::GetNumGenes(void)
{
  return m_NumGenes;
}/* end GetNumGenes() */

/******************************************************************************
SetMutationRate()

Sets the mutation rate of all genes in the Gene List.
******************************************************************************/
void Chromosome::SetMutationRate(double rate)
{
   int i;

   for(i = 0; i < m_NumGenes; i++){ m_pGenes[i]->SetMutationRate(rate); }
} /* end SetMutationRate() */

/******************************************************************************
Crossover()

Performs crossover between the genes of this chromosome class and the genes of
the pMate argument. The crossed-over genes will replace the genes of this 
chromosome.
******************************************************************************/
void Chromosome::Crossover(Chromosome * pMate)
{
   Gene * pMom;
   Gene * pPop;
   int i;

   double F1 = this->GetFitness();
   double F2 = pMate->GetFitness();

   for(i = 0; i < m_NumGenes; i++)
   {
      pMom = m_pGenes[i];
      pPop = pMate->GetGenePtr(i);
      pMom->Crossover(pPop, F1, F2, m_NumGenes);
   }/* end for() */
}/* end Crossover() */

/******************************************************************************
Mutate()

Mutates the genes of this chromosome according to the pre-established mutation 
rate. Every time a mutation occurs, the 'pCount' array is updated, so that
mutation metrics can be tracked.
******************************************************************************/
void Chromosome::Mutate(int * pCount)
{
   int i;
   for(i = 0; i < m_NumGenes; i++){ pCount[i] += m_pGenes[i]->Mutate(); }
}/* end Mutate() */

/******************************************************************************
Destroy()

Frees up memory used by the Genes of this chromosome.
******************************************************************************/
void Chromosome::Destroy(void)
{
   int i;
  
   for(i = 0; i < m_NumGenes; i++) 
   { 
      delete m_pGenes[i];
   }
   delete [] m_pGenes;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Copy()

Copies this chromosome data in 'pCopy' into the member variables of the 'this'
class. If the number of genes in 'pCopy' is not the same as the number of genes
in 'this', then no action is taken.
******************************************************************************/
void Chromosome::Copy(Chromosome * pCopy)
{
   int i;
   Gene * pGene;      
   
   if(m_NumGenes != pCopy->GetNumGenes()){ return;}

   m_Fitness = pCopy->GetFitness();

   for(i = 0; i < m_NumGenes; i++)
   {
      pGene = pCopy->GetGenePtr(i);
      m_pGenes[i]->Copy(pGene);
   }/* end for() */
}/* end Copy() */

/******************************************************************************
CreateRandomChromo()

Creates a chromosome whose Genes have randomly assigned values. Returns a 
pointer to this chromosome.
******************************************************************************/
Chromosome *  Chromosome::CreateRandomChromo(void)
{
   Chromosome * pChromo;
   Gene * pRandom;
   int i;

   NEW_PRINT("Chromosome", 1);
   pChromo = new Chromosome(m_Fitness, m_NumGenes);   
   MEM_CHECK(pChromo);

   for(i = 0; i < m_NumGenes; i++)
   {
      pRandom = m_pGenes[i]->CreateRandomGene();
      pChromo->SetGenePtr(pRandom, i);      
   }/* end for() */

   return pChromo;
}/* end CreateRandomChromo() */

/******************************************************************************
CreateChromo()

Creates a chromosome whose Genes have values assigned based on the vals array.

Returns a pointer to this chromosome.
******************************************************************************/
Chromosome *  Chromosome::CreateChromo(double * vals)
{
   Chromosome * pChromo;
   Gene * pGene;
   int i;

   NEW_PRINT("Chromosome", 1);
   pChromo = new Chromosome(m_Fitness, m_NumGenes);   
   MEM_CHECK(pChromo);

   for(i = 0; i < m_NumGenes; i++)
   {
      pGene = m_pGenes[i]->CreateGene(vals[i]);
      pChromo->SetGenePtr(pGene, i);
   }/* end for() */

   return pChromo;
}/* end CreateChromo() */

