/******************************************************************************
File     : Gene.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

A Gene is an encoded design variable. A sequence of Genes is the major compnent
of the Chromosome, which in turn makes up the contents of a ChromosomePool. 
Various Genetic Algorithm operations can be performed on a Gene, including 
Random Instantiation, Crossover, Mutation and Cloning.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes, metrics collection
10-19-05    lsm   Replaced rand() with MyRand()
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Gene.h"

#include "Utility.h"
#include "Exception.h"

/******************************************************************************
CTOR

Constructs a Gene using the real number arg. and its upper and lower bounds 
and the mutation rate.
******************************************************************************/
RealEncodedGene::RealEncodedGene(double val, double lwr, double upr){
   m_Value        = val;
   m_LowerBound   = lwr;
   m_UpperBound   = upr;


   IncCtorCount();
} /* end RealEncodedGene::CTOR */


/******************************************************************************
CTOR

Constructs an empty Gene
******************************************************************************/
RealEncodedGene::RealEncodedGene(void) {
    m_Value = NULL;
    m_LowerBound = NULL;
    m_UpperBound = NULL;

    IncCtorCount();
} /* end RealEncodedGene::CTOR */



/******************************************************************************
GetRandomValue()

Generates a random value between lower and upper bound of the gene
******************************************************************************/
double RealEncodedGene::GetRandomValue(void) {
    //generate a random between lower and upper bound
    double r, val, range;
    range = m_UpperBound - m_LowerBound;
    r = (double)MyRand() / (double)MY_RAND_MAX;
    val = (r * range) + m_LowerBound;

    return val;
} 


/******************************************************************************
Copy()

Creates a copy of the 'pCopy' gene.
******************************************************************************/
void RealEncodedGene::Copy(Gene * pCopy)
{
   m_Value = pCopy->GetValue();
   m_LowerBound = pCopy->GetLwr();
   m_UpperBound = pCopy->GetUpr();

} /* end RealEncodedGene::Copy() */

/******************************************************************************
CreateRandomGene()

Generates a random gene between lower and upper bound. 

Returns a pointer to the gene.
******************************************************************************/
Gene * RealEncodedGene::CreateRandomGene(void)
{
   Gene * pGene;
   //generate a random between lower and upper bound
   double r, val, range;
   range = m_UpperBound - m_LowerBound;
   r = (double)MyRand() / (double)MY_RAND_MAX;
   val = (r * range) + m_LowerBound;

   pGene = new RealEncodedGene(val, m_LowerBound, m_UpperBound);   
   MEM_CHECK(pGene);

   return pGene;
} /* end RealEncodedGene::CreateRandomGene() */

/******************************************************************************
CreateGene()

Generates a gene using the given value.

Returns a pointer to the gene.
******************************************************************************/
Gene * RealEncodedGene::CreateGene(double val)
{
   Gene * pGene; 

   NEW_PRINT("RealEncodedGene", 1);
   pGene = new RealEncodedGene(val, m_LowerBound, m_UpperBound);   
   MEM_CHECK(pGene);

   return pGene;
} /* end RealEncodedGene::CreateGene() */

