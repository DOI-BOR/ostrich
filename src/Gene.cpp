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
RealEncodedGene::RealEncodedGene(double val, double lwr, double upr,  double rate, double xover){
   m_Value        = val;
   m_LowerBound   = lwr;
   m_UpperBound   = upr;
   m_MutationRate = rate;
   m_CrossoverRate = xover;

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
    m_MutationRate = NULL;
    m_CrossoverRate = NULL;

    IncCtorCount();
} /* end RealEncodedGene::CTOR */

/******************************************************************************
Crossover()

Performs crossover using the convex arithmetic crossover technique. The child 
is a weigthed average of the parents, such that the child value always lies 
between those of its parents. The weights used are based on the fitness values of
the parents. A random epsilon value between +/- 10% is also added.
******************************************************************************/
double RealEncodedGene::Crossover(RealEncodedGene* pMate, double selfValue, double mateVal, double F1, double F2, int np)
{
   double r, w1, w2, p, s;
   double range, lwr, upr;
   double childVal;

   //determine range of gene
   lwr = this->GetLwr(); upr = this->GetUpr(); range = upr -lwr;

   //determine weights, based on fitness values
   p = 1.00-(MyMin(fabs(F1),fabs(F2))/MyMax(fabs(F1),fabs(F2)));
   if(CheckOverflow(p)) p = 0.00;

   if(F1 > F2) {
      w1=MyMin(0.5+0.5*p,1.00);
      w2=1.00-w1;
   }
   else {
      w2=MyMin(0.5+0.5*p,1.00);
      w1=1.00-w2;
   }

   r = (double)MyRand() / (double)MY_RAND_MAX;

   if (r < m_CrossoverRate) {      
      childVal = (selfValue * w1) + (mateVal * w2);

      r = (double)MyRand() / (double)MY_RAND_MAX; //0 to 1
      s = (double)MyRand() / (double)MY_RAND_MAX; //0 to 1
      r *= 0.20; //0 to 20%
      r -= 0.10; //-10% to +10%

      // epsilon perturbation using normal distribution
      // centered on childVal with standard deviation estimated 
      // by the fitness value
      double sd = sqrt(fabs(MyMax(F1,F2))/(double)np);
      childVal = MyGaussRand(childVal, sd);

/*
      if(s < 0.0) //value-relative epsilon
      {
         childVal *= (1.00+r);
      }
      else //range-relative epsilon
      {
         childVal += (range*r*0.1);
      }
*/

      //enforce parameter limits
      if(childVal > upr) childVal = selfValue + (upr - selfValue)*s;
      if(childVal < lwr) childVal = selfValue - (selfValue - lwr)*s;
        
      return childVal;
   } 
   else {
       return selfValue;
   }
}/* end RealEncodedGene::Crossover() */

/******************************************************************************
Mutate()

Mutates the gene at random. If mutation occurs, the gene is assigned a random
value between the upper and lower bound.

Returns 0 if no mutation occurs, a 1 if a mutation did occur.
******************************************************************************/
double RealEncodedGene::Mutate(double inputValue)
{
   double r;
   double range;

   range = m_UpperBound - m_LowerBound;

   r = (double)MyRand() / (double)MY_RAND_MAX;

   if(r < m_MutationRate)
   {      
      r = (double)MyRand() / (double)MY_RAND_MAX;
      m_Value = (r * range) + m_LowerBound;
      return m_Value;
   }
   else {
       return inputValue;
   }
   
} 


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
   m_MutationRate = pCopy->GetMutationRate();
   m_CrossoverRate = pCopy->GetCrossoverRate();
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

   pGene = new RealEncodedGene(val, m_LowerBound, m_UpperBound, m_MutationRate, m_CrossoverRate);   
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
   pGene = new RealEncodedGene(val, m_LowerBound, m_UpperBound, m_MutationRate, m_CrossoverRate);   
   MEM_CHECK(pGene);

   return pGene;
} /* end RealEncodedGene::CreateGene() */

