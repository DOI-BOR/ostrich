/******************************************************************************
File     : Chromosome.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

A Chromosome is a set of design variables (also called genes) that make up a 
single solution to a given optimization problem. As their name implies, 
chromosomes are used by the Genetic Algorithm.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
07-08-04    lsm   moved DTOR to Destroy() so CTOR/DTOR bookeeping is correct
08-17-04    lsm   RAM fragmentation fixes
******************************************************************************/
#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "MyHeaderInc.h"

//forward declarations
class Gene;

/******************************************************************************
class Chromosome
******************************************************************************/
class Chromosome
{
   private:      
      double m_Fitness;
      Gene ** m_pGenes;
      int m_NumGenes;

   public:
      Chromosome(double fitness, int numGenes);
      ~Chromosome(void){ DBG_PRINT("Chromosome::DTOR"); Destroy();}
      void Destroy(void);

      void Crossover(Chromosome * pMate);
      void Mutate(int * pCount);

      int  GetNumGenes(void);

      Gene * GetGenePtr(int i);
      void SetGenePtr(Gene * pGene, int i);

      void SetFitness(double fitness);
      double GetFitness(void);

      void Copy(Chromosome * pCopy);
      Chromosome * CreateRandomChromo(void);
      Chromosome * CreateChromo(double * vals);

      void SetMutationRate(double rate);      
}; /* end class Chromosome */

#endif /* CHROMOSOME_H */

