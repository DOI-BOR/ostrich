/******************************************************************************
File     : BinaryGene.cpp
Author   : L. Shawn Matott
Copyright: 2005, L. Shawn Matott

A Gene is an encoded design variable. A sequence of Genes is the major compnent
of the Chromosome, which in turn makes up the contents of a ChromosomePool. 
Various Genetic Algorithm operations can be performed on a Gene, including 
Random Instantiation, Crossover, Mutation and Cloning.

Version History
10-18-05    lsm   created
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Gene.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
CTOR

Constructs a Gene using the number arg. and its upper and lower bounds 
and the mutation rate.
******************************************************************************/
BinaryEncodedGene::BinaryEncodedGene(double val, double lwr, double upr, 
                                 double rate, double xover)
{
   m_Range = (int)(upr-lwr);
   m_Offset = (int)lwr;
   m_NumBits = CalcNumBits(m_Range);
   if(m_NumBits > 32) 
   {
      LogError(ERR_BGA, "BinaryEncodedGene(): coding exceeds 32 bits");
      ExitProgram(1);
   }   
   m_BitMask = CalcBitMask(m_NumBits);
   m_Value = ((int)(val-m_Offset))&(m_BitMask);
   m_LowerBound   = lwr;
   m_UpperBound   = upr;
   m_MutationRate = rate;
   m_CrossoverRate = xover;

   IncCtorCount();
} /* end BinaryEncodedGene::CTOR */

/******************************************************************************
Crossover()

Performs crossover using the bit-flipping technique. All bits following a 
randomly generated bit position are interchagned. The crossover value replaces 
the value of the parent (i.e. 'this').
******************************************************************************/
void BinaryEncodedGene::Crossover(Gene * pMate, double F1, double F2, int np)
{
   int r; //randomly generated bit position
   int mask; //bit mask
   int childVal;
   int mateVal;
   double rd;

   rd = (double)MyRand() / (double)MY_RAND_MAX;

   if(rd < m_CrossoverRate)
   {      
      mateVal = ((BinaryEncodedGene *)pMate)->GetCodedValue();
      r = (MyRand() % m_NumBits)+1;
      mask = CalcBitMask(r);

      mateVal  = (mateVal & mask);
      childVal = (m_Value & ~mask);
      m_Value = childVal | mateVal;
      if(m_Value > m_Range) m_Value = m_Range;
   }
}/* end BinaryEncodedGene::Crossover() */

/******************************************************************************
Mutate()

Mutates the gene at random. If mutation occurs, the gene is assigned a random
value between the upper and lower bound.

Returns 0 if no mutation occurs, a 1 if a mutation did occur.
******************************************************************************/
int BinaryEncodedGene::Mutate(void)
{
   double p; //randomly generated mutation probability
   int r; //randomly generated bit position
   int mask; //bit mask
   int childVal;
   int mateVal;
  
   p = (double)MyRand() / (double)MY_RAND_MAX;

   if(p < m_MutationRate)
   {      
      mateVal = (MyRand() % (m_Range+1));
      r = (MyRand() % m_NumBits)+1;
      mask = CalcBitMask(r);

      mateVal  = (mateVal & mask);
      childVal = (m_Value & ~mask);
      m_Value = childVal | mateVal;
      if(m_Value > m_Range) m_Value = m_Range;

      return 1;
   }/* end if() */
   return 0;
} /* end BinaryEncodedGene::Mutate() */

/******************************************************************************
Copy()

Creates a copy of the 'pCopy' gene.
******************************************************************************/
void BinaryEncodedGene::Copy(Gene * pCopy)
{
   m_Value      = ((BinaryEncodedGene *)pCopy)->GetCodedValue();
   m_LowerBound = pCopy->GetLwr();
   m_UpperBound = pCopy->GetUpr();
   m_MutationRate = pCopy->GetMutationRate();
   m_CrossoverRate = pCopy->GetCrossoverRate();
} /* end BinaryEncodedGene::Copy() */

/******************************************************************************
CreateRandomGene()

Generates a random gene between lower and upper bound. 

Returns a pointer to the gene.
******************************************************************************/
Gene * BinaryEncodedGene::CreateRandomGene(void)
{
   Gene * pGene;
   //generate a random between lower and upper bound
   double r, val, range;
   range = m_UpperBound - m_LowerBound;
   r = (double)MyRand() / (double)MY_RAND_MAX;
   val = (r * range) + m_LowerBound;

   NEW_PRINT("BinaryEncodedGene", 1);
   pGene = new BinaryEncodedGene(val, m_LowerBound, m_UpperBound, m_MutationRate, m_CrossoverRate);   
   MEM_CHECK(pGene);

   return pGene;
} /* end BinaryEncodedGene::CreateRandomGene() */

/******************************************************************************
CreateGene()

Generates a gene using the given value.

Returns a pointer to the gene.
******************************************************************************/
Gene * BinaryEncodedGene::CreateGene(double val)
{
   Gene * pGene; 

   NEW_PRINT("BinaryEncodedGene", 1);
   pGene = new BinaryEncodedGene(val, m_LowerBound, m_UpperBound, m_MutationRate, m_CrossoverRate);   
   MEM_CHECK(pGene);

   return pGene;
} /* end BinaryEncodedGene::CreateGene() */

/******************************************************************************
CalcNumBits()

Compute the number of bits required to represent the given range.
******************************************************************************/
int BinaryEncodedGene::CalcNumBits(int range)
{
   int i;

   for(i = 0; i <= 32; i++)
   {
      if((int)(pow((double)2,i)) > range) break;
   }

   return i;   
}/* end CalcNumBits() */

/******************************************************************************
CalcBitMask()

Calculate a mask that will clear all bits after 'pos' bits. If pos = 8 then 
mask will be: 0000 0000 0000 0000 0000 0000 1111 1111 bin = 0x000000FF hex

******************************************************************************/
int BinaryEncodedGene::CalcBitMask(int pos)
{
   int i, mask = 0x00000000;

   for(i = 0; i < pos; i++)
   {
      mask |= 1 << i;
   }

   return mask;
}/* end CalcBitMask() */
