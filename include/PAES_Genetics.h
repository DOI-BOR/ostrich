/******************************************************************************
File       : PAES_Genetics.h
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)

PAES - Pareto archived evolution strategy

Version History
12-29-17    lsm   created file
******************************************************************************/
#ifndef PAES_GENETICS_H
#define PAES_GENETICS_H

#include "PAES_Helpers.h"
#include "PAES_ObjFunc.h"

/******************************************************************************
class PAES_Gene
 
Abstract class that represents a gene usable by the PAES algorithm
******************************************************************************/
class PAES_Gene
{
   public:
      PAES_Random * m_random; // Random number generator
      PAES_VariableType m_geneType; // Type of the gene
      FILE * m_cout;

      // Constructors
      PAES_Gene(PAES_VariableType geneType, PAES_Random * random, FILE * pFile);
      PAES_Gene(PAES_Gene & gene);
      PAES_Gene(PAES_Gene * gene);
  
      // Destructor
      virtual ~PAES_Gene(void);
  
      // Methods
      virtual int bitFlipMutation(double mutationProbability); 
  
      virtual int randomMutation(double mutationProbability);
      virtual int polynomialMutation(double mutationProbability, double distributionIndex); 
      virtual int uniformMutation(double mutationProbability, double perturbation); 

      virtual double getRealAllele(void);
      virtual void   writeGenotype(FILE * pFile) = 0;
  
      void print(void); 
  
      // Operators 
      PAES_Gene & operator=(const PAES_Gene& gene);
}; /* end class PAES_Gene */

/******************************************************************************
 class PAES_BinaryGene
 ******************************************************************************/
class PAES_BinaryGene : public PAES_Gene 
{
   public:
      char * m_allele; // Bit string
      int    m_numberOfBits; // Number of bits of the bit string
  
      // Constructors  
      PAES_BinaryGene(double X, int bits, PAES_Random * random, FILE * pFile);
      PAES_BinaryGene(int bits, PAES_Random * random, FILE * pFile);
      PAES_BinaryGene(PAES_BinaryGene & binaryGene);
      PAES_BinaryGene(PAES_BinaryGene * binaryGene);

      // Destructors
      ~PAES_BinaryGene(void);
  
      // Methods
      int bitFlipMutation(double mutationProbability);

      void writeGenotype(FILE * pFile);
  
      void print(void);

      // Operators
      PAES_BinaryGene & operator=(const PAES_BinaryGene& binaryGene);
}; /* end class PAES_BinaryGene */

/******************************************************************************
class PAES_BinaryRealGene
******************************************************************************/
class PAES_BinaryRealGene : public PAES_Gene 
{
   private:
      double decodeToReal(void);

   public:
      char * m_binaryAllele; // Binary string
      int    m_numberOfBits; // Number of bits of the bit string
      double m_realAllele; // Real value of the allele
      double m_lowerBound; // Lower bound of the allele
      double m_upperBound; // Upper bound of the allele
  
      // Constructors
      PAES_BinaryRealGene(double X, int numberOfBits, double lowerBound, double upperBound, PAES_Random * random, FILE * pFile);
      PAES_BinaryRealGene(int numberOfBits, double lowerBound, double upperBound, PAES_Random * random, FILE * pFile);
      PAES_BinaryRealGene(PAES_BinaryRealGene & binaryRealGene);
      PAES_BinaryRealGene(PAES_BinaryRealGene * binaryRealGene);

      // Destructor
      virtual ~PAES_BinaryRealGene(void);
  
      // Methods
      int    bitFlipMutation(double mutationProbability);
      double getRealAllele(void);
      void   writeGenotype(FILE * pFile);
      void   encodeFromReal(double X);

      void print(void);

      // Operators
      PAES_BinaryRealGene & operator=(const PAES_BinaryRealGene& gene);
}; /* end class PAES_BinaryRealGene */

/******************************************************************************
class PAES_RealGene
******************************************************************************/
class PAES_RealGene : public PAES_Gene 
{
   public:
      double m_allele; // Allele
      double m_lowerBound; // Lower bound of the allele
      double m_upperBound; // Upper bound of the allele

      // Constructors
      PAES_RealGene(PAES_Random * random, FILE * pFile);
      PAES_RealGene(double X, double lowerBound, double upperBound, PAES_Random * random, FILE * pFile);
      PAES_RealGene(double lowerBound, double upperBound, PAES_Random * random, FILE * pFile);
      PAES_RealGene(PAES_RealGene & realgene);
      PAES_RealGene(PAES_RealGene * realGene);

      // Destructor
      ~PAES_RealGene(void);
  
      // Methods
      int randomMutation(double mutationProbability); 
      int polynomialMutation(double mutationProbability, double distributionIndex); 
      int uniformMutation(double mutationProbability, double perturbation); 

      double getRealAllele(void);
      void   writeGenotype(FILE * pFile);

      void print(void);

      // Operators
      PAES_RealGene & operator=(const PAES_RealGene& realGene);
}; /* end class PAES_RealGene */

/******************************************************************************
class PAES_Chromosome
******************************************************************************/
class PAES_Chromosome 
{
   public:
      int m_length; // Chromosome lenght
      PAES_MultiobjectiveProblem * m_problem; // Problem to solve 
      PAES_Gene ** m_gene; // Genes of the chromosome
      FILE * m_cout;

      // Constructors
      PAES_Chromosome(double * X, PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile); 
      PAES_Chromosome(PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile);
      PAES_Chromosome(PAES_Chromosome & chromosome);
      PAES_Chromosome(PAES_Chromosome * chromosome);
  
      // Destructor
      ~PAES_Chromosome(void);
  
      void print(void);

      // Operators
      PAES_Chromosome & operator=(PAES_Chromosome & chromosome);
}; /* end class Chromosome */

/******************************************************************************
class PAES_Individual
******************************************************************************/
class PAES_Individual 
{
   public:
      PAES_MultiobjectiveProblem * m_problem; // Problem to be solved
      PAES_Random * m_random; // Random number generator
      PAES_Chromosome * m_chromosome; // Chromosome of the individual
      double * m_fitness; // Array of fitness values
      FILE * m_cout;
  
      int m_gridLocation; // Necessary if an adaptive grid is used
  
      // Constructors
      PAES_Individual(PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile);
      PAES_Individual(PAES_Individual & individual);
      PAES_Individual(PAES_Individual * individual);
      PAES_Individual(double * X, double * F, PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile);
  
      // Destructor
      ~PAES_Individual(void);
  
      // Methods
      void    setFitness(double * fitness);
      double *getFitness(void) const;  
  
      int  dominanceTest(PAES_Individual * individual);
      int  constraintComparison(PAES_Individual * individual);
      bool identicalFitness(PAES_Individual * individual);
    
      // Binary mutation operators
      int bitFlipMutation(double mutationProbability);
  
      // Real mutation operators
      int randomMutation(double mutationProbability);
      int polynomialMutation(double mutationProbability, double distributionIndex);
      int uniformMutation(double mutationProbability, double perturbation);

      void print(void);
  
      // Operators
      PAES_Individual & operator=(PAES_Individual &individual);
  
      void printFitness(FILE * pOut);
}; /* end class PAES_Individual */

/******************************************************************************
class PAES_Population
 
This class represents a population of individuals. The population size can be d
ynamic: individuals can be deleted, and new individuals can be added. 
The number of elements is bounded.
******************************************************************************/
class PAES_Population 
{
   protected:
      int           m_populationSize; // The number of individuals
      int           m_maximumPopulationSize; // The maximun number of individuals
      PAES_Individual ** m_population; // The vector of individuals
      PAES_Random      * m_random; // Random number management
      PAES_MultiobjectiveProblem * m_problem; // Problem to solve
      FILE * m_cout;

   public:
      // Constructor
      PAES_Population(int populationSize, int maximumPopulationSize, PAES_Random * random,
                      PAES_MultiobjectiveProblem * problem, FILE * pFile);
             
      // Destructor
      ~PAES_Population(void);

      // Methods
      int getPopulationSize(void) const;
      int getMaximumPopulationSize(void) const;
      PAES_Individual * getIth(int index) const;
      void setIth(int index, PAES_Individual * ind);
      void deleteIth(int index);
      void addIndividual(PAES_Individual * individual);
      void setFitness(int index, double * fitness);
  
      void printFitness(char  * fileName);
      void printGenotype(char * fileName);
}; /* end class PAES_Population */

/******************************************************************************
class PAES_AdaptiveGrid
******************************************************************************/
class PAES_AdaptiveGrid 
{
   public:
      int    * m_hypercube; // Hypercube division for keeping diversity
      double * m_divisionSize; // Division sizes of the adaptive grid
      double * m_gridLimits; // Limits of the adaptive grid
      long     m_currentGridSize; // Current size of the adaptive grid
      int      m_mostCrowdedHypercube;

      int m_numberOfFunctions; // 
      int m_depth; // 
    
      double * m_upperBestFitness;
      double * m_lowerBestFitness;
  
      // Constructor
      PAES_AdaptiveGrid(void);
      PAES_AdaptiveGrid(int depth, int numberOfFunctions);
  
      // Destructor
      ~PAES_AdaptiveGrid(void);
  
      // Methods
      void updateGridLocations(PAES_Population * population, PAES_Individual * individual);
      int  findLocation(PAES_Individual * individual);
  
   private:
      // These variables are used in findLocation(). They are defined here for 
      // efficiency purposes
      int    * m_increment;
      double * m_tmpDivisionSize;
}; /* end class PAES_AdaptiveGrid */

#endif /* PAES_GENETICS_H */
