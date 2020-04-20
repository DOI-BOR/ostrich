/******************************************************************************
File     : Gene.h
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
10-19-05    lsm   Added support for BGA
******************************************************************************/
#ifndef GENE_H
#define GENE_H

#include "MyHeaderInc.h"

/******************************************************************************
class Gene

 Represents a charcteristic/parameter of a design.

 The class is an abstract class.

 The following methods define the characteristics of this class:
	Crossover(), Mutate() and GetValue() 
******************************************************************************/
class Gene
{   
   public:
      virtual ~Gene(void){ DBG_PRINT("Gene::DTOR"); }
      virtual void Destroy(void)=0;

      virtual void Crossover(Gene * pMate, double F1, double F2, int np)=0;
      virtual int  Mutate   (void)        =0;

      virtual double GetValue(void)       = 0;
      virtual void   SetValue(double val) = 0;
            
      virtual void   Copy            (Gene * pCopy) = 0;
      virtual Gene * CreateRandomGene(void)         = 0;
      virtual Gene * CreateGene      (double val)   = 0;

      virtual double GetUpr(void) = 0;
      virtual double GetLwr(void) = 0;

      virtual double GetMutationRate(void) = 0;
      virtual void   SetMutationRate(double rate) = 0;

      virtual double GetCrossoverRate(void) = 0;
      virtual void   SetCrossoverRate(double rate) = 0;
}; /* end class Gene */

/******************************************************************************
class RealEncodedGene

Represents a  gene encoded in the real number form. Here the gene is 
reprsented as a real number.

The variables lowerBound and upperBound are values between which
the value of the gene is allowed to fluctuate.

The upper and lower bound are used so that we have control over mutation. 
This allows us to ensure that the mutation does not change its value
to an infeasible one.
******************************************************************************/
class  RealEncodedGene : public Gene
{   
   public:
      RealEncodedGene(double val, double lwr, double upr, double rate, double xover);
     ~RealEncodedGene(void){ DBG_PRINT("RealEncodedGene::DTOR"); Destroy(); }

     void Destroy(void) { IncDtorCount();}
     
     void Crossover(Gene * pMate, double F1, double F2, int np);
     int Mutate    (void);

     void Copy              (Gene * pCopy);
     Gene * CreateRandomGene(void);
     Gene * CreateGene      (double val);
  
     void   SetValue(double val) { m_Value = val;}
     double GetValue(void)       { return m_Value;}    
     double GetUpr  (void)       { return m_UpperBound;}
     double GetLwr  (void)       { return m_LowerBound;}

     double GetMutationRate(void)       { return m_MutationRate;}
     void   SetMutationRate(double rate){ m_MutationRate = rate;}

     double GetCrossoverRate(void)       { return m_CrossoverRate;}
     void   SetCrossoverRate(double rate){ m_CrossoverRate = rate;}

   private:      
      double m_Value;
      double m_LowerBound;
      double m_UpperBound;

   protected:                                    
      double m_MutationRate;
      double m_CrossoverRate;

}; /* end class RealEncodedGene */

/******************************************************************************
class BinaryEncodedGene

Represents a  gene encoded in binary form. Here the gene is reprsented as a 
string of bits.

The variables lowerBound and upperBound are values between which
the value of the gene is allowed to fluctuate.

The upper and lower bound are used so that we have control over mutation. 
This allows us to ensure that the mutation does not change its value
to an infeasible one.
******************************************************************************/
class  BinaryEncodedGene : public Gene
{   
   public:
      BinaryEncodedGene(double val, double lwr, double upr, double rate, double xover);
     ~BinaryEncodedGene(void){ DBG_PRINT("BinaryEncodedGene::DTOR"); Destroy(); }

     void Destroy(void) { IncDtorCount();}
     
     void Crossover(Gene * pMate, double F1, double F2, int np);
     int Mutate    (void);

     void Copy              (Gene * pCopy);
     Gene * CreateRandomGene(void);
     Gene * CreateGene      (double val);
  
     void   SetValue(double val) { m_Value = ((int)(val-m_Offset))&(m_BitMask);}
     double GetValue(void)       { return (double)(m_Value+m_Offset);}
     double GetUpr  (void)       { return m_UpperBound;}
     double GetLwr  (void)       { return m_LowerBound;}
     double GetMutationRate(void)       { return m_MutationRate;}
     void   SetMutationRate(double rate){ m_MutationRate = rate;}
     double GetCrossoverRate(void)       { return m_CrossoverRate;}
     void   SetCrossoverRate(double rate){ m_CrossoverRate = rate;}

   private:      
      int GetCodedValue(void) { return m_Value;}
      int CalcNumBits(int range);
      int CalcBitMask(int pos);
      double m_LowerBound;
      double m_UpperBound;
      int m_Value;   //bit representation of values
      int m_Offset;  //lower bound, cast as an int
      int m_Range;   //upper-lower
      int m_NumBits; //bits required to represent full range
      int m_BitMask; //mask that zeros unused bits

   protected:                                    
      double m_MutationRate;
      double m_CrossoverRate;
}; /* end class BinaryEncodedGene */

#endif /* GENE_H */


