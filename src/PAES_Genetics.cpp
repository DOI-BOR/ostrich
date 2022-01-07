/******************************************************************************
File       : PAES_Genetics.cpp
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)

PAES - Pareto archived evolution strategy

Version History
12-29-17    lsm   created file
******************************************************************************/
#include "PAES_Genetics.h"
 
/******************************************************************************
PAES_Gene::CTOR
******************************************************************************/
PAES_Gene::PAES_Gene(PAES_VariableType geneType, PAES_Random * random, FILE * pFile)
{
   m_random = random;
   m_geneType = geneType;
   m_cout = pFile;
} /* end PAES_Gene::CTOR */

/******************************************************************************
PAES_Gene::CTOR

Copy constructor.
******************************************************************************/
PAES_Gene::PAES_Gene(PAES_Gene & gene) 
{
   m_random = gene.m_random;
   m_geneType = gene.m_geneType;
   m_cout = gene.m_cout;
} /* end PAES_Gene copy constructor */

/******************************************************************************
PAES_Gene::CTOR

Copy constructor.
******************************************************************************/
PAES_Gene::PAES_Gene(PAES_Gene * gene) 
{
   m_random = gene->m_random;
   m_geneType = gene->m_geneType;
   m_cout = gene->m_cout;
} /* end PAES_Gene copy constructor */

/******************************************************************************
PAES_Gene::DTOR
******************************************************************************/
PAES_Gene::~PAES_Gene(void) 
{
} /* end PAES_Gene::DTOR */

/******************************************************************************
PAES_Gene::bitFlipMutation()
******************************************************************************/
int PAES_Gene::bitFlipMutation(double mutationProbability) 
{
   fprintf(m_cout, "Bit-flip mutation cannot be applied to a gene of type %d \n", m_geneType); 
   fclose(m_cout);
   exit(-1);
} /* and PAES_Gene::bitFlipMutation() */

/******************************************************************************
PAES_Gene::randomMutation()
******************************************************************************/
int PAES_Gene::randomMutation(double mutationProbability) 
{
   fprintf(m_cout, "Random mutation cannot be applied to a gene of type %d\n", m_geneType); 
   fclose(m_cout); 
   exit(-1);
} /* end PAES_Gene::randomMutation() */

/******************************************************************************
PAES_Gene::polynomialMutation()
******************************************************************************/
int PAES_Gene::polynomialMutation(double mutationProbability, double distributionIndex) 
{
   fprintf(m_cout, "Polynomial mutation cannot be applied to a gene of type %d\n", m_geneType); 
   fclose(m_cout); 
   exit(-1);
} /* end PAES_Gene::polynomialMutation() */ 

/******************************************************************************
PAES_Gene::uniformMutation()
******************************************************************************/
int PAES_Gene::uniformMutation(double mutationProbability, double perturbation) 
{
   fprintf(m_cout, "Uniform mutation cannot be applied to a gene of type %d\n", m_geneType); 
   fclose(m_cout); 
   exit(-1);
} /* end PAES_Gene::uniformMutation() */

/******************************************************************************
PAES_Gene::getRealAllele()
******************************************************************************/
double PAES_Gene::getRealAllele(void) 
{
   fprintf(m_cout, "getRealAllele() cannot be applied to a gene of type %d\n", m_geneType); 
   fclose(m_cout);
   exit(-1);
} /* end PAES_Gene::getRealAllele() */

/******************************************************************************
PAES_Gene::operator=
******************************************************************************/
PAES_Gene & PAES_Gene::operator=(const PAES_Gene& gene) 
{
   m_random = gene.m_random;
   m_geneType = gene.m_geneType;  
   return *this;  
} /* end Gene::operator= */

/******************************************************************************
PAES_Gene::print
******************************************************************************/
void PAES_Gene::print(void) 
{
   fprintf(m_cout, "Gene type: %d\n", m_geneType);
} /* end PAES_Gene::print() */ 

/******************************************************************************
PAES_BinaryGene::CTOR
******************************************************************************/
PAES_BinaryGene::PAES_BinaryGene(double X, int numberOfBits, PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_BINARY, random, pFile)
{
   int i, pow_bit, iX;
  
   m_numberOfBits = numberOfBits;
   m_cout = pFile;
   m_allele = new char[numberOfBits];
  
   if (m_allele == NULL)
   {
      fprintf(m_cout, "BinaryGene::BinaryGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   iX = (int)(X + 0.5);
   for (i = 0; i < numberOfBits; i++)
   {
      pow_bit = (int)pow(2, i);
      m_allele[i] = '0' + (iX & pow_bit) / pow_bit;
   }
} /* end PAES_BinaryGene::CTOR */

/******************************************************************************
PAES_BinaryGene::CTOR
******************************************************************************/
PAES_BinaryGene::PAES_BinaryGene(int numberOfBits, PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_BINARY, random, pFile)
{
   int i;
  
   m_numberOfBits = numberOfBits;
   m_cout = pFile;
   m_allele = new char[numberOfBits];
  
   if (m_allele == NULL)
   {
      fprintf(m_cout, "BinaryGene::BinaryGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   for (i = 0; i < numberOfBits; i++)
   {
      if (m_random->rnd(0,1) == 1)
      {
         m_allele[i] = '1';
      }
      else
      {
         m_allele[i] = '0';
      }
   }
} /* end PAES_BinaryGene::CTOR */

/******************************************************************************
PAES_BinaryGene::Copy CTOR
******************************************************************************/
PAES_BinaryGene::PAES_BinaryGene(PAES_BinaryGene & binaryGene) : PAES_Gene(binaryGene) 
{
   int i;
  
   m_numberOfBits= binaryGene.m_numberOfBits;
   m_allele = new char[m_numberOfBits];

   if (m_allele == NULL) 
   {
      fprintf(m_cout, "BinaryGene::BinaryGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end  if() */

   for (i = 0; i < m_numberOfBits; i++)
   {
      m_allele[i] = binaryGene.m_allele[i];
   }  
} /* end PAES_BinaryGene::copy CTOR */

/******************************************************************************
PAES_BinaryGene::copy CTOR
******************************************************************************/
PAES_BinaryGene::PAES_BinaryGene(PAES_BinaryGene * binaryGene) : PAES_Gene(binaryGene) 
{
   int i;
  
   m_numberOfBits = binaryGene->m_numberOfBits;
   m_allele = new char[m_numberOfBits];

   if (m_allele == NULL) 
   {
      fprintf(m_cout, "BinaryGene::BinaryGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */

   for (i = 0; i < m_numberOfBits; i++)
   {
      m_allele[i] = binaryGene->m_allele[i];
   }  
} /* end PAES_BinaryGene::copy CTOR */

/******************************************************************************
PAES_BinaryGene::DTOR
******************************************************************************/
PAES_BinaryGene::~PAES_BinaryGene(void) 
{
   delete [] m_allele;
} /* end PAES_BinaryGene::DTOR */

/******************************************************************************
PAES_BinaryGene::bitFlipMutation()
******************************************************************************/
int PAES_BinaryGene::bitFlipMutation(double mutationProbability) 
{
   int mutations;
  
   mutations = 0;

   for (int i = 0; i < m_numberOfBits; i++)
   {
      if (m_random->flip((float)mutationProbability) == 1) 
      {
         mutations++;
         if (m_allele[i] == '1')
         {
            m_allele[i] = '0';
         }
         else
         {
            m_allele[i] = '1';
         }
      } /* end if() */
   }
   return mutations;
} /* end PAES_BinaryGene::bitFlipMutation() */

/******************************************************************************
PAES_BinaryGene::writeGenotype()
******************************************************************************/
void PAES_BinaryGene::writeGenotype(FILE * pFile) 
{
   int i;
   for (i = 0; i < m_numberOfBits; i++)
   {
      fprintf(pFile, "%c", m_allele[i]);
   }/* end for() */
} /* end PAES_BinaryGene::writeGenotype() */

/******************************************************************************
PAES_BinaryGene::operator=
******************************************************************************/
PAES_BinaryGene & PAES_BinaryGene::operator=(const PAES_BinaryGene& binaryGene) 
{
   int i;

   m_numberOfBits = binaryGene.m_numberOfBits;
   for (i = 0; i < m_numberOfBits; i++)
   {
      m_allele[i] = binaryGene.m_allele[i];  
   }  
   return *this;
} /* end PAES_BinaryGene::operator= */

/******************************************************************************
PAES_BinaryGene::print()
******************************************************************************/
void PAES_BinaryGene::print(void) 
{
   int i;
  
   PAES_Gene::print();

   fprintf(m_cout, " Bits: %d\n", m_numberOfBits);
   fprintf(m_cout, " allele: \n");
   for (i = 0; i < m_numberOfBits; i++)
   {
      if (m_allele[i] == '1')
      {  
         fprintf(m_cout, "1");
      }
      else
      {
         fprintf(m_cout, "0");  
      }
   }/* end for() */
} /* end PAES_BinaryGene::print() */ 

/******************************************************************************
PAES_BinaryRealGene::CTOR
******************************************************************************/
PAES_BinaryRealGene::PAES_BinaryRealGene(double X, int numberOfBits, double lowerBound, double upperBound, PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_BINARY_REAL, random, pFile) 
{                    
   m_upperBound = upperBound;
   m_lowerBound = lowerBound;                    
   m_numberOfBits = numberOfBits;
   m_cout = pFile;
  
   m_binaryAllele = new char[numberOfBits];
  
   if (m_binaryAllele == NULL) 
   {
      fprintf(m_cout, "BinaryRealGene::BinaryRealGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   encodeFromReal(X);
} /* end PAES_BinaryRealGene::BinaryRealGene() */

/******************************************************************************
PAES_BinaryRealGene::CTOR
******************************************************************************/
PAES_BinaryRealGene::PAES_BinaryRealGene(int numberOfBits, double lowerBound, double upperBound, PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_BINARY_REAL, random, pFile) 
{
   int i;                    
                    
   m_upperBound = upperBound;
   m_lowerBound = lowerBound;                    
   m_numberOfBits = numberOfBits;
   m_cout = pFile;
  
   m_binaryAllele = new char[numberOfBits];
  
   if (m_binaryAllele == NULL) 
   {
      fprintf(m_cout, "BinaryRealGene::BinaryRealGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   for (i = 0; i < numberOfBits; i++)
   {
      if (m_random->rnd(0,1) == 1)
      {
         m_binaryAllele[i] = '1';
      }
      else
      {
         m_binaryAllele[i] = '0';
      }
   }/* end for() */

   m_realAllele = decodeToReal();      
} /* end PAES_BinaryRealGene::BinaryRealGene() */

/******************************************************************************
PAES_BinaryRealGene::copy CTOR
******************************************************************************/
PAES_BinaryRealGene::PAES_BinaryRealGene(PAES_BinaryRealGene & gene) : PAES_Gene(gene) 
{
   int i;                    

   m_upperBound = gene.m_upperBound;
   m_lowerBound = gene.m_lowerBound;                                                      
   m_numberOfBits = gene.m_numberOfBits;
   m_binaryAllele = new char[m_numberOfBits];
  
   if (m_binaryAllele == NULL) 
   {
      fprintf(m_cout, "BinaryRealGene::BinaryRealGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   for (i = 0; i < m_numberOfBits; i++)
   {
      m_binaryAllele[i] = gene.m_binaryAllele[i];
   }
    
   m_realAllele = gene.m_realAllele;
} /* end PAES_BinaryRealGene::BinaryRealGene() */

/******************************************************************************
PAES_BinaryRealGene::copy CTOR
******************************************************************************/
PAES_BinaryRealGene::PAES_BinaryRealGene(PAES_BinaryRealGene * gene) : PAES_Gene(gene) 
{
   int i;                    
                    
   m_upperBound = gene->m_upperBound;
   m_lowerBound = gene->m_lowerBound;                                                      
   m_numberOfBits = gene->m_numberOfBits;
   m_cout = gene->m_cout;
   m_binaryAllele = new char[m_numberOfBits];
  
   if (m_binaryAllele == NULL) 
   {
      fprintf(m_cout, "BinaryRealGene::BinaryRealGene-> Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   for (i = 0; i < m_numberOfBits; i++)
   {
      m_binaryAllele[i] = gene->m_binaryAllele[i];
   }
    
   m_realAllele = gene->m_realAllele;
} /* end PAES_BinaryRealGene::BinaryRealGene() */

/******************************************************************************
PAES_BinaryRealGene::DTOR
******************************************************************************/
PAES_BinaryRealGene::~PAES_BinaryRealGene(void) 
{
   delete [] m_binaryAllele;
} /* end PAES_BinaryRealGene::~BinaryRealGene() */

/******************************************************************************
PAES_BinaryRealGene::operator=
******************************************************************************/
PAES_BinaryRealGene & PAES_BinaryRealGene::operator=(const PAES_BinaryRealGene&  gene) 
{
   int i;                    

   m_upperBound = gene.m_upperBound;
   m_lowerBound = gene.m_lowerBound;                          
   m_numberOfBits = gene.m_numberOfBits;
   for (i = 0; i < m_numberOfBits; i++)
   {
      m_binaryAllele[i] = gene.m_binaryAllele[i];
   }
   m_realAllele = gene.m_realAllele;        
  
   return *this;                                             
} /* end PAES_BinaryRealGene::operator= */                                                    

/******************************************************************************
PAES_BinaryRealGene::bitFlipMutation()
******************************************************************************/
int PAES_BinaryRealGene::bitFlipMutation(double mutationProbability) 
{
   int mutations;
   int i;
   mutations = 0;

   for (i = 0; i < m_numberOfBits; i++)
   {
      if (m_random->flip((float)mutationProbability) == 1) 
      {
         mutations ++;
         if (m_binaryAllele[i] == '1')
         {
            m_binaryAllele[i] = '0';
         }
         else
         {
            m_binaryAllele[i] = '1';
         }
      } /* end if() */
   }/* end for() */
    
   m_realAllele = decodeToReal();      
   return mutations;
} /* end PAES_BinaryRealGene::bitFlipMutation() */

/******************************************************************************
PAES_BinaryRealGene::getRealAllele()
******************************************************************************/
double PAES_BinaryRealGene::getRealAllele(void) 
{
   return m_realAllele;
} /* end PAES_BinaryRealGene::getRealAllele() */

/******************************************************************************
PAES_BinaryRealGene::writeGenotype()
******************************************************************************/
void PAES_BinaryRealGene::writeGenotype(FILE * pFile) 
{
   fprintf(pFile, "%E", m_realAllele);
} /* end PAES_BinaryRealGene::writeGenotype() */

/******************************************************************************
PAES_BinaryRealGene::print()
******************************************************************************/
void PAES_BinaryRealGene::print(void) 
{
   int i;
   PAES_Gene::print();

   fprintf(m_cout, "Real allele: %E", m_realAllele);
  
   fprintf(m_cout, " Bits: %d\n", m_numberOfBits);
   fprintf(m_cout, " Binary allele: ");
   for (i = 0; i < m_numberOfBits; i++)
   {
      fprintf(m_cout, "%c", m_binaryAllele[i]);
   }/* end for() */
   fprintf(m_cout, "\n");
} /* end PAES_BinaryRealGene::print() */ 

/******************************************************************************
PAES_BinaryRealGene::decodeToReal()
******************************************************************************/
double PAES_BinaryRealGene::decodeToReal(void) 
{
   int    i;
   double realValue;
  
   realValue = 0.0;
   for (i = 0; i < m_numberOfBits; i++) 
   {
      if (m_binaryAllele[i] == '1')
      {
         realValue += pow((double)2.0, i);
      }
   }/* end for() */
    
   return (m_lowerBound + realValue*(m_upperBound - m_lowerBound) /(pow(2.0, m_numberOfBits) - 1));
} /* end PAES_BinaryRealGene::decodeToReal() */

/******************************************************************************
PAES_BinaryRealGene::encodeFromReal()
******************************************************************************/
void PAES_BinaryRealGene::encodeFromReal(double X) 
{
   int    i;
   double frac, range;
   int pow_frac, denom, pow_bit;

   denom = (int)(pow(2, m_numberOfBits) - 1);
   range = m_upperBound - m_lowerBound;

   frac = (X - m_lowerBound) / range;
   pow_frac =  (int)(0.5 + (frac * denom));
  
   for (i = 0; i < m_numberOfBits; i++) 
   {
      pow_bit = (int)pow(2, i);
      m_binaryAllele[i] = '0' + (pow_frac & pow_bit) / pow_bit;
   }/* end for() */

   m_realAllele = X;
} /* end PAES_BinaryRealGene::encodeFromReal() */

/******************************************************************************
PAES_RealGene::CTOR
******************************************************************************/
PAES_RealGene::PAES_RealGene(PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_REAL, random, pFile) 
{
   m_lowerBound = PAES_MIN_REAL;
   m_upperBound = PAES_MAX_REAL;
   m_cout = pFile;
   m_allele = m_random->rndreal((float)m_lowerBound, (float)m_upperBound);
} /* end PAES_RealGene::PAES_RealGene() */

/******************************************************************************
PAES_RealGene::CTOR
******************************************************************************/
PAES_RealGene::PAES_RealGene(double lowerBound, double upperBound, PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_REAL, random, pFile) 
{
   m_lowerBound = lowerBound;
   m_upperBound = upperBound;
   m_allele = m_random->rndreal((float)lowerBound, (float)upperBound);
} /* end PAES_RealGene::PAES_RealGene() */

/******************************************************************************
PAES_RealGene::CTOR
******************************************************************************/
PAES_RealGene::PAES_RealGene(double X, double lowerBound, double upperBound, PAES_Random * random, FILE * pFile) : PAES_Gene(PAES_REAL, random, pFile) 
{
   m_lowerBound = lowerBound;
   m_upperBound = upperBound;
   m_allele = X;
} /* end PAES_RealGene::PAES_RealGene() */

/******************************************************************************
PAES_RealGene::copy CTOR
******************************************************************************/
PAES_RealGene::PAES_RealGene(PAES_RealGene & realGene) : PAES_Gene(realGene) 
{
   m_allele     = realGene.m_allele;
   m_lowerBound = realGene.m_lowerBound;
   m_upperBound = realGene.m_upperBound;  
} /* end PAES_RealGene::PAES_RealGene() */

/******************************************************************************
PAES_RealGene::copy CTOR
******************************************************************************/
PAES_RealGene::PAES_RealGene(PAES_RealGene * realGene) : PAES_Gene(realGene) 
{
   m_allele     = realGene->m_allele;
   m_lowerBound = realGene->m_lowerBound;
   m_upperBound = realGene->m_upperBound;
} /* end PAES_RealGene::PAES_RealGene() */

/******************************************************************************
PAES_RealGene::DTOR
******************************************************************************/
PAES_RealGene::~PAES_RealGene(void) 
{
} /* end PAES_RealGene::~PAES_RealGene() */

/******************************************************************************
PAES_RealGene::randomMutation()
******************************************************************************/
int PAES_RealGene::randomMutation(double mutationProbability) 
{
   int    mutations;
   double rnd;
  
   mutations = 0;
  
   rnd = m_random->rndreal(0.0, 1.0);
   if (rnd <= mutationProbability) 
   {
      m_allele = m_lowerBound + rnd *(m_upperBound - m_lowerBound);
      mutations++;
   } /* end if() */    
  
   return mutations;
} /* end PAES_RealGene::uniformMutation() */

/******************************************************************************
PAES_RealGene::polynomialMutation()
******************************************************************************/
int PAES_RealGene::polynomialMutation(double mutationProbability, double distributionIndex) 
{                                   
   double temp;
   int    mutations;
   double rnd;

   double delta;
   double deltaq;
   double mu;
  
   mutations = 0;

   rnd = m_random->rndreal(0.0, 1.0);

   if (rnd <= mutationProbability) 
   {
      // calculate delta
      if (m_allele > m_lowerBound) 
      { 
         if ((m_allele - m_lowerBound) < (m_upperBound - m_allele))
         {
            delta = (m_allele - m_lowerBound) / (m_upperBound - m_allele);
         }
         else
         {
            delta = (m_upperBound - m_allele) / (m_upperBound - m_lowerBound);
         }
         rnd = m_random->rndreal(0.0, 1.0);
         mu  = 1.0/(distributionIndex + 1.0);
         if(rnd <= 0.5) 
         {
            double xy = 1.0-delta;
            temp = 2*rnd+(1-2*rnd)*(pow(xy,(distributionIndex + 1)));
            deltaq =  pow(temp, mu) - 1.0;
         } /* end if() */
         else 
         {
            double xy = 1.0-delta;
            temp = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(distributionIndex + 1)));
            deltaq = 1.0 - (pow(temp,mu));
         } /* end else() */

         /*Change the value of the allele */
         m_allele += deltaq * (m_upperBound - m_lowerBound);
         if (m_allele < m_lowerBound)
         {
            m_allele = m_lowerBound;
         }
         if (m_allele > m_upperBound)
         {
            m_allele = m_upperBound;
         }
      } /* end if() */
      else 
      {
         rnd = m_random->rndreal(0.0, 1.0);
         m_allele = rnd * (m_upperBound - m_lowerBound) + m_lowerBound;
      } /* end else() */
      mutations ++;
   } /* end if() */

   return mutations;
} /* end PAES_RealGene::polynomialMutation() */

/******************************************************************************
PAES_RealGene::uniformMutation()
******************************************************************************/
int PAES_RealGene::uniformMutation(double mutationProbability, double perturbation) 
{
   int    mutations;
   double rnd;
  
   mutations = 0;
  
   rnd = m_random->rndreal(0.0, 1.0);
   if (rnd <= mutationProbability) 
   {
      double tmp = (rnd + 0.5)*perturbation;
      m_allele = m_allele + tmp;
      if (m_allele < m_lowerBound)
      {
         m_allele = m_lowerBound;
      }
      if (m_allele > m_upperBound)
      {
         m_allele = m_upperBound;
      }
   
      mutations ++;
   } /* end if() */
  
   return mutations;
} /* end PAES_RealGene::uniformMutation() */

/******************************************************************************
PAES_RealGene::getRealAllele()
******************************************************************************/
double PAES_RealGene::getRealAllele(void) 
{
   return m_allele;
} /* end PAES_RealGene::getRealAllele() */

/******************************************************************************
PAES_RealGene::writeGenotype()
******************************************************************************/
void PAES_RealGene::writeGenotype(FILE * pFile) 
{
   fprintf(pFile, "%E", m_allele);
} /* end PAES_RealGene::writeGenotype() */

/******************************************************************************
PAES_RealGene::operator=
******************************************************************************/
PAES_RealGene & PAES_RealGene::operator=(const PAES_RealGene& realGene) 
{
   m_allele = realGene.m_allele;

   return *this;  
} /* end PAES_RealGene::operator= */

/******************************************************************************
PAES_RealGene::print()
******************************************************************************/
void PAES_RealGene::print(void) 
{
   PAES_Gene::print();
   fprintf(m_cout, " allele: %E\n", m_allele);
} /* end PAES_RealGene::print() */ 

/******************************************************************************
PAES_Chromosome::CTOR
******************************************************************************/
PAES_Chromosome::PAES_Chromosome(double * X, PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile) 
{
   m_problem = problem;
   m_length  = problem->m_numberOfVariables;
   m_cout = pFile;

   m_gene = new PAES_Gene*[m_length];
  
   if (m_gene == NULL) 
   {
      fprintf(m_cout, "Chromosome::Chromosome->Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   int i;
   for (i = 0; i < m_length; i++) 
   {
      switch (m_problem->m_variable[i]) 
      {
         case PAES_BINARY:
            m_gene[i] = new PAES_BinaryGene(X[i], problem->m_numberOfBits, random, pFile);
            break;

         case PAES_REAL:
            m_gene[i] = new PAES_RealGene(X[i], m_problem->m_lowerLimit[i], m_problem->m_upperLimit[i], random, m_cout);
            break;

         case PAES_BINARY_REAL:
            m_gene[i] = new PAES_BinaryRealGene(X[i], problem->m_bitsPerVariable[i], m_problem->m_lowerLimit[i], m_problem->m_upperLimit[i], random, m_cout);
            break;

         default: 
            fprintf(m_cout, "Chromosome::Chromosome->variable type %d unknown\n", m_problem->m_variable[i]);
            fclose(m_cout);
            exit(-1);
      } /* end switch() */
   } /* end for() */
} /* end PAES_Chromosome::PAES_Chromosome() */

/******************************************************************************
PAES_Chromosome::CTOR
******************************************************************************/
PAES_Chromosome::PAES_Chromosome(PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile) 
{
   m_problem = problem;
   m_length  = problem->m_numberOfVariables;
   m_cout = pFile;

   m_gene = new PAES_Gene*[m_length];
  
   if (m_gene == NULL) 
   {
      fprintf(m_cout, "Chromosome::Chromosome->Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   int i;
   for (i = 0; i < m_length; i++) 
   {
      switch (m_problem->m_variable[i]) 
      {
         case PAES_BINARY:
            m_gene[i] = new PAES_BinaryGene(problem->m_numberOfBits, random, pFile);
            break;

         case PAES_REAL:
            m_gene[i] = new PAES_RealGene(m_problem->m_lowerLimit[i], m_problem->m_upperLimit[i], random, m_cout);
            break;

         case PAES_BINARY_REAL:
            m_gene[i] = new PAES_BinaryRealGene(problem->m_bitsPerVariable[i], m_problem->m_lowerLimit[i], m_problem->m_upperLimit[i], random, m_cout);
            break;

         default: 
            fprintf(m_cout, "Chromosome::Chromosome->variable type %d unknown\n", m_problem->m_variable[i]);
            fclose(m_cout);
            exit(-1);
      } /* end switch() */
   } /* end for() */
} /* end PAES_Chromosome::PAES_Chromosome() */

/******************************************************************************
PAES_Chromosome::copy CTOR
******************************************************************************/
PAES_Chromosome::PAES_Chromosome(PAES_Chromosome & chromosome) 
{
   m_problem = chromosome.m_problem;
   m_length  = chromosome.m_length;
   m_gene = new PAES_Gene*[m_length];

   if (m_gene == NULL) 
   {
      fprintf(m_cout, "Chromosome::Chromosome->Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */

   int i;
   for (i = 0; i < m_length; i++) 
   {
      switch (m_problem->m_variable[i]) 
      {
         case PAES_BINARY:
            m_gene[i] = new PAES_BinaryGene((PAES_BinaryGene*)chromosome.m_gene[i]);
            break;

         case PAES_REAL:
            m_gene[i] = new PAES_RealGene((PAES_RealGene*)chromosome.m_gene[i]);
            break;

         case PAES_BINARY_REAL:
            m_gene[i] = new PAES_BinaryRealGene((PAES_BinaryRealGene*)chromosome.m_gene[i]);
            break;

         default: 
            fprintf(m_cout, "Chromosome::Chromosome->variable type %d unknown\n", m_problem->m_variable[i]);
            fclose(m_cout);
            exit(-1);
      } /* end switch */
   } /* end for() */
} /* end PAES_Chromosome::PAES_Chromosome() */

/******************************************************************************
PAES_Chromosome::copy CTOR
******************************************************************************/
PAES_Chromosome::PAES_Chromosome(PAES_Chromosome * chromosome) 
{
   m_problem = chromosome->m_problem;
   m_length  = chromosome->m_length;
   m_gene = new PAES_Gene*[m_length];

   if (m_gene == NULL) 
   {
      fprintf(m_cout, "Chromosome::Chromosome->Error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */

   int i;
   for (i = 0; i < m_length; i++) 
   {
      switch (m_problem->m_variable[i]) 
      {
         case PAES_BINARY:
            m_gene[i] = new PAES_BinaryGene((PAES_BinaryGene*)chromosome->m_gene[i]);
            break;

         case PAES_REAL:
            m_gene[i] = new PAES_RealGene((PAES_RealGene*)chromosome->m_gene[i]);
            break;

         case PAES_BINARY_REAL:
            m_gene[i] = new PAES_BinaryRealGene((PAES_BinaryRealGene*)chromosome->m_gene[i]);
            break;

         default: 
            fprintf(m_cout, "Chromosome::Chromosome->variable type %d unknown\n", m_problem->m_variable[i]);
            fclose(m_cout);
            exit(-1);
      } /* end switch */
   } /* end for() */
} /* end PAES_Chromosome::PAES_Chromosome() */

/******************************************************************************
PAES_Chromosome::DTOR
******************************************************************************/
PAES_Chromosome::~PAES_Chromosome(void) 
{
   int i;

   for (i = 0; i < m_length; i++)
   {
      switch (m_problem->m_variable[i]) 
      {
         case PAES_BINARY:
            delete (PAES_BinaryGene *)m_gene[i];
            break;

         case PAES_REAL:
            delete (PAES_RealGene *)m_gene[i];
            break;

         case PAES_BINARY_REAL:
            delete (PAES_BinaryRealGene *)m_gene[i];
            break;

         default: 
            fprintf(m_cout, "Chromosome::~Chromosome->variable type %d unknown\n", m_problem->m_variable[i]);
            fclose(m_cout);
            exit(-1);
      } /* end switch */
   }/* end for() */
    
   delete [] m_gene;  
} /* end PAES_Chromosome::~PAES_Chromosome() */

/******************************************************************************
PAES_Chromosome::operator=
******************************************************************************/
PAES_Chromosome & PAES_Chromosome::operator=(PAES_Chromosome & chromosome) 
{
   m_problem = chromosome.m_problem;
   m_length  = chromosome.m_length;
   int i;

   for (i = 0; i < m_length; i++) 
   {
      switch (m_problem->m_variable[i]) 
      {
         case PAES_BINARY:
            *((PAES_BinaryGene *)m_gene[i]) = *((PAES_BinaryGene*)chromosome.m_gene[i]);
            break;

         case PAES_REAL:
            *((PAES_RealGene *)m_gene[i]) = *((PAES_RealGene*)chromosome.m_gene[i]);
            break;

         case PAES_BINARY_REAL:
            *((PAES_BinaryRealGene *)m_gene[i]) = *((PAES_BinaryRealGene*)chromosome.m_gene[i]);
            break;

         default: 
            fprintf(m_cout, "Chromosome->operator= -> variable type %d unknown\n", m_problem->m_variable[i]);
            fclose(m_cout);
            exit(-1);
      } /* end switch */
   } /* end for() */
   return chromosome;
} /* end PAES_Chromosome::operator= */

/******************************************************************************
PAES_Chromosome::print()
******************************************************************************/
void PAES_Chromosome::print(void) 
{
   int i;
  
   fprintf(m_cout, "Length: %d\n", m_length);

   for (i = 0; i < m_length; i++)
   {
      if (m_gene[i]->m_geneType == PAES_BINARY)
      {
         ((PAES_BinaryGene*)(m_gene[i]))->print();
      }
      else if (m_gene[i]->m_geneType == PAES_REAL)
      {
         ((PAES_RealGene*)(m_gene[i]))->print();
      }
      else if (m_gene[i]->m_geneType == PAES_BINARY_REAL)
      {
         ((PAES_BinaryRealGene*)(m_gene[i]))->print();
      }
      else 
      { 
         fprintf(m_cout, "Chromosome: variable type %d unknown\n", m_gene[i]->m_geneType);
         exit(-1);
      } /* end else() */
   }/* end for() */
   fprintf(m_cout, "\n");
} /* end PAES_Chromosome::print() */ 

/******************************************************************************
PAES_Individual::CTOR
******************************************************************************/
PAES_Individual::PAES_Individual(PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile) 
{
   int i;
  
   m_problem    = problem; 
   m_random     = random;
   m_cout = pFile;
   m_chromosome = new PAES_Chromosome(m_problem, m_random, m_cout);
   m_fitness    = new double[m_problem->m_numberOfFunctions];
  
   if ((m_fitness == NULL) || (m_chromosome == NULL)) 
   {
      fprintf(m_cout, "Individual::Individual() error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   for (i = 0; i <m_problem->m_numberOfFunctions; i++)
   {
      m_fitness[i] = 0.0;
   }
} /* end PAES_Individual::PAES_Individual() */

/******************************************************************************
PAES_Individual::CTOR
******************************************************************************/
PAES_Individual::PAES_Individual(double * X, double * F, PAES_MultiobjectiveProblem * problem, PAES_Random * random, FILE * pFile) 
{
   int i;
  
   m_problem    = problem; 
   m_random     = random;
   m_cout = pFile;
   m_chromosome = new PAES_Chromosome(X, m_problem, m_random, m_cout);
   m_fitness    = new double[m_problem->m_numberOfFunctions];
  
   if ((m_fitness == NULL) || (m_chromosome == NULL)) 
   {
      fprintf(m_cout, "Individual::Individual() error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   for (i = 0; i < m_problem->m_numberOfFunctions; i++)
   {
      m_fitness[i] = F[i];
   }
} /* end PAES_Individual::PAES_Individual() */

/******************************************************************************
PAES_Individual::copy CTOR
******************************************************************************/
PAES_Individual::PAES_Individual(PAES_Individual & individual) 
{
   int i;
  
   m_problem    = individual.m_problem;
   m_chromosome = new PAES_Chromosome(individual.m_chromosome);
   m_fitness    = new double[individual.m_problem->m_numberOfFunctions];
  
   if ((m_fitness == NULL) || (m_chromosome == NULL)) 
   {
      fprintf(m_cout, "Individual::Individual() error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */

   for (i = 0; i <m_problem->m_numberOfFunctions; i++)
   {
      m_fitness[i] = individual.m_fitness[i];
   }
} /* end PAES_Individual::PAES_Individual() */

/******************************************************************************
PAES_Individual::copy CTOR
******************************************************************************/
PAES_Individual::PAES_Individual(PAES_Individual * individual) 
{
   int i;
  
   m_problem    = individual->m_problem;
   m_chromosome = new PAES_Chromosome(individual->m_chromosome);
   m_fitness    = new double[individual->m_problem->m_numberOfFunctions];
  
   if ((m_fitness == NULL) || (m_chromosome == NULL)) 
   {
      fprintf(m_cout, "Individual::Individual() -> error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */

   for (i = 0; i <m_problem->m_numberOfFunctions; i++)
   {
      m_fitness[i] = individual->m_fitness[i];
   }
} /* end PAES_Individual::PAES_Individual() */

/******************************************************************************
PAES_Individual::DTOR
******************************************************************************/
PAES_Individual::~PAES_Individual(void) 
{
   delete m_chromosome;  
   delete [] m_fitness;
} /* end PAES_Individual::~PAES_Individual() */

/******************************************************************************
PAES_Individual::setFitness()

Sets the fitness of the individual

fitness: An array of doubles with the fitness of the individual
******************************************************************************/
void PAES_Individual::setFitness(double * fitness) 
{
   int i;
   for (i = 0; i < m_problem->m_numberOfFunctions; i++) 
   {
      m_fitness[i] = fitness[i];
   } /* end for() */
} /* end PAES_Individual::setFitness() */

/******************************************************************************
PAES_Individual::getFitness()

Gets the fitness of the individual

Returns an array of doubles with the fitness of the individual
******************************************************************************/
double * PAES_Individual::getFitness(void) const 
{
   return m_fitness;
} /* end PAES_Individual::getFitness() */

/******************************************************************************
PAES_Individual::bitFlipMutation()
******************************************************************************/
int PAES_Individual::bitFlipMutation(double mutationProbability) 
{
   int i;
   int mutations; 
  
   mutations = 0;
   for (i = 0; i < m_chromosome->m_length; i++) 
   {
      mutations += m_chromosome->m_gene[i]->bitFlipMutation(mutationProbability);
   } /* end for() */
  
   return mutations;
} /* end PAES_Individual::bitFlipMutation() */

/******************************************************************************
PAES_Individual::uniformMutation()
******************************************************************************/
int PAES_Individual::uniformMutation(double mutationProbability, double perturbation) 
{
   int i;
   int mutations; 
  
   mutations = 0;
   for (i = 0; i < m_chromosome->m_length; i++) 
   {
      mutations += m_chromosome->m_gene[i]->uniformMutation(mutationProbability, perturbation);
   } /* end for() */
  
   return mutations;
} /* end PAES_Individual::uniformMutation() */

/******************************************************************************
PAES_Individual::randomMutation()
******************************************************************************/
int PAES_Individual::randomMutation(double mutationProbability) 
{
   int i;
   int mutations; 
  
   mutations = 0;
   for (i = 0; i < m_chromosome->m_length; i++) 
   {
      mutations += m_chromosome->m_gene[i]->randomMutation(mutationProbability);
   } /* end for() */
  
   return mutations;
} /* end PAES_Individual::randomMutation() */

/******************************************************************************
PAES_Individual::polynomialMutation()

Applies a polynomial mutation with certain probability.

individual: The individual to mutate

mutationProbability: The probability of a bit mutation

Returns The number of mutations
 
 The polynomial mutation is defined in [Deb 2002], p 124.
******************************************************************************/
int PAES_Individual::polynomialMutation(double mutationProbability, double distributionIndex) 
{
   int i;
   int mutations; 
   
   mutations = 0;
   for (i = 0; i < m_chromosome->m_length; i++) 
   {
      mutations += m_chromosome->m_gene[i]->polynomialMutation(mutationProbability, distributionIndex);
   } /* end for() */
  
   return mutations;                                   
} /* end PAES_Individual::polynomialMutation() */

/******************************************************************************
PAES_Individual::identicalFitness()
******************************************************************************/
bool PAES_Individual::identicalFitness(PAES_Individual * individual) 
{
   int i;

   for (i = 0; i < m_problem->m_numberOfFunctions; i++)
   {
      if (m_fitness[i] != individual->m_fitness[i])
      {
         return false;
      }
   }

   return true;  
} /* end PAES_Individual::identicalFitness() */

/******************************************************************************
PAES_Individual::constraintComparison()

Constraint comparison test between two individuals.

individual: The individual against the test if performed

Returns: 1 if this individual violates less constraints than the individual 
           passed as parameter
        -1 if this individual violates less constraints than the current 
           individual
         0 otherwise
It is assumed that the individual has been previously decoded to a real
representation.
******************************************************************************/
int PAES_Individual::constraintComparison(PAES_Individual * individual) 
{
   int constraintsFirstIndividual;
   int constraintsSecondIndividual;        
  
   int result;
  
   result = 0;
   constraintsFirstIndividual  = m_problem->numberOfNonSatisfiedConstraints(this);
   constraintsSecondIndividual = m_problem->numberOfNonSatisfiedConstraints(individual);
   if (constraintsFirstIndividual < constraintsSecondIndividual)
   {
      result = 1;   // Current dominates 
   }
   else if (constraintsFirstIndividual > constraintsSecondIndividual)
   {
      result = -1;  // Current is dominated
   }
   return result;                                     
} /* end PAES_Individual::constraintComparison() */                                  

/******************************************************************************
PAES_Individual::dominanceTest()

individual: The individual against which the test is performed

Returns: +1 if the individual passed as argument is dominated
         -1 if this individual dominates the current individual
          0 if the two individuals are non-dominated
******************************************************************************/
int PAES_Individual::dominanceTest(PAES_Individual * individual) 
{
   int i;                        
   int last;                        
   int current; 
   int result; 
  
   bool finished;

   finished = false;
 
   if (m_problem->m_numberOfConstraints > 0) 
   {
      result = this->constraintComparison(individual);
      if (result != 0)
      {
         finished = true;
      }
   } /* end if() */
 
   last = 0;
   i    = 0; 
   
   while (!finished) 
   {
      if (this->m_fitness[i] < individual->m_fitness[i])
      {
         current = 1;
      }
      else if (this->m_fitness[i] > individual->m_fitness[i])
      {
         current = -1;
      }
      else
      {
         current = 0;
      }
    
      if ((current != 0) && (current == -last)) 
      {
         finished = true;
         result   = 0;
      } /* end if() */
      else 
      {
         last = current;
      
         i ++;
         if (i == this->m_problem->m_numberOfFunctions) 
         {
            finished = true;    
            result   = last;
         } /* end if() */
      } /* end else() */
   } /* end while() */
  
   return result;  
} /* end PAES_Individual::dominanceTest() */

/******************************************************************************
PAES_Individual::operator=
******************************************************************************/
PAES_Individual & PAES_Individual::operator=(PAES_Individual &individual) {
   int i;
  
   m_problem     = individual.m_problem;
   *m_chromosome = *(individual.m_chromosome);
  
   for (i = 0; i <m_problem->m_numberOfFunctions; i++)
   {
      m_fitness[i] = individual.m_fitness[i];  
   }

   return individual;
} /* end PAES_Individual::operator= */

/******************************************************************************
PAES_Individual::print()
******************************************************************************/
void PAES_Individual::print(void) 
{
   int i;
   m_chromosome->print();

   fprintf(m_cout, "Fitness: ");
   for (i = 0; i < m_problem->m_numberOfFunctions; i++)
   {
      fprintf(m_cout, "%E ", m_fitness[i]); 
   }
   fprintf(m_cout, "\n");
} /* end PAES_Individual::print() */ 

/******************************************************************************
PAES_Individual::printFitness()
******************************************************************************/
void PAES_Individual::printFitness(FILE * pOut) 
{
   int i;
   fprintf(pOut, "Fitness: ");
   for (i = 0; i < m_problem->m_numberOfFunctions; i++)
   {
      fprintf(pOut, "%E ", m_fitness[i]);
   }
   fprintf(pOut, "\n");
} /* end PAES_Individual::printFitness() */

/******************************************************************************
PAES_Population::CTOR

populationSize: The size of the population
chromosomeLength: The length of the individual chromosomes
numberOfFunctions: The number of objective functions
******************************************************************************/
PAES_Population::PAES_Population(int populationSize, int maximumPopulationSize, PAES_Random * random, PAES_MultiobjectiveProblem * problem, FILE * pFile) 
{
   m_problem               = problem;
   m_populationSize        = populationSize;
   m_maximumPopulationSize = maximumPopulationSize;
   m_random                = random;
   m_cout = pFile;
  
   m_population = new PAES_Individual *[maximumPopulationSize];
   if (!m_population) 
   {
      fprintf(m_cout, "Population::Population() -> error when asking for memory\n");
      fclose(m_cout);
      exit(-1);
   } /* end if() */
  
   int i;
   for (i = 0; i < m_populationSize; i++) 
   {
      m_population[i] = new PAES_Individual(m_problem, random, m_cout);
   }
} /* end PAES_Population::Population() */

/******************************************************************************
PAES_Population::DTOR
******************************************************************************/
PAES_Population::~PAES_Population(void) 
{
   return;
} /* end PAES_Population::~Population() */

/******************************************************************************
PAES_Population::getPopulationSize()
******************************************************************************/
int PAES_Population::getPopulationSize(void) const 
{
   return m_populationSize;
} /* end PAES_Population::getPopulationSize() */
  
/******************************************************************************
PAES_Population::getMaximumPopulationSize()
******************************************************************************/
int PAES_Population::getMaximumPopulationSize(void) const 
{
   return m_maximumPopulationSize;
} /* end PAES_Population::getMaximumPopulationSize() */

/******************************************************************************
PAES_Population::getIth()

Gets the i-th individual of the population

index: The index of the individual
******************************************************************************/
PAES_Individual * PAES_Population::getIth(int index) const 
{
   if ((index < m_populationSize) && (index >= 0))
   {
      return m_population[index];
   }
   else 
   {
      fprintf(m_cout, "Population::getIth - Index out of range when getting a copy of an individual\n");
      fclose(m_cout);
      exit(-1);
   } /* end else() */
} /* end PAES_Population::getIth() */

/******************************************************************************
PAES_Population::setIth()

Sets the i-th individual of the population

index: The index of the individual
individual: The individual to assign
******************************************************************************/
void PAES_Population::setIth(int index, PAES_Individual * individual) 
{
   if ((index < m_populationSize) && (index >= 0))
   {
      m_population[index] = individual;
   }
   else 
   {
      fprintf(m_cout, "Population::setIth() - Index out of range when inserting an individual\n");
      fclose(m_cout);
      exit(-1);
   } /* end else() */
} /* end PAES_Population::setIth() */

/******************************************************************************
PAES_Population::setFitness()

Sets the fitness of the i-th individual

index: The position of the individual
fitness: The individual fitness
******************************************************************************/
void PAES_Population::setFitness(int index, double * fitness) 
{
   m_population[index]->setFitness(fitness);
} /* end PAES_Population::setFitness() */

/******************************************************************************
PAES_Population::addIndividual()
******************************************************************************/
void PAES_Population::addIndividual(PAES_Individual * individual) 
{
   if (m_populationSize < m_maximumPopulationSize) 
   {
      m_population[m_populationSize] = individual;
      m_populationSize += 1;
   } /* end if() */
   else 
   {
      fprintf(m_cout, "Population::addIndividual() -> The population size has reach to its maximum size\n");
      fclose(m_cout);
      exit(-1);
   } /* end else() */
} /* end PAES_Population::addIndividual() */

/******************************************************************************
PAES_Population::deleteIth()
******************************************************************************/
void PAES_Population::deleteIth(int index) 
{
   if ((index < m_populationSize) && (index >= 0)) 
   {
      delete m_population[index];
      int i;
      for (i = index; i < (m_populationSize - 1); i++)
      {
         m_population[i] = m_population[i + 1];
      }
      m_populationSize--;
   } /* end if() */
   else 
   {
      fprintf(m_cout, "Population::deleteIth() - Index %d out of range\n", index);
      fclose(m_cout);
      exit(-1);
   } /* end else() */
} /* end PAES_Population::deleteIth() */

/******************************************************************************
PAES_Population::printFitness()

Prints the fitness of the individuals in the population

fileName: The name of the ouput file
******************************************************************************/
void PAES_Population::printFitness(char *fileName) 
{
   FILE * outputFile;
   int      i;
   int      j;

   outputFile = fopen(fileName, "w");
  
   for (j = 0; j < getPopulationSize(); j++) 
   {
      for (i = 0; i < m_problem->m_numberOfFunctions; i++)
      {
         fprintf(outputFile, "%E\t", getIth(j)->getFitness()[i]);;
      }
      fprintf(outputFile, "\n");;
   } /* end for() */

   fclose(outputFile);      
} /* end PAES_Population::printFitness() */

/******************************************************************************
PAES_Population::printGenotype()

Prints the genotype of the individuals in the population

fileName: The name of the ouput file
******************************************************************************/
void PAES_Population::printGenotype(char *fileName) 
{
   FILE * outputFile;
   int      i;
   int      j;

   outputFile = fopen(fileName, "w");
  
   for (j = 0; j < getPopulationSize(); j++) 
   {
      for (i = 0; i < m_problem->m_numberOfVariables; i++) 
      {
         getIth(j)->m_chromosome->m_gene[i]->writeGenotype(outputFile);
         fprintf(outputFile, "\t");
      } /* end for() */
      fprintf(outputFile, "\n");
   } /* end for() */
  
   fclose(outputFile);      
} /* end PAES_Population::printGenotype() */

/******************************************************************************
PAES_AdaptiveGrid::CTOR
******************************************************************************/
PAES_AdaptiveGrid::PAES_AdaptiveGrid(void) 
{
   m_numberOfFunctions = 0;
   m_depth = 0;
   m_currentGridSize = 0;
   m_hypercube = NULL;
   m_divisionSize = NULL;
   m_gridLimits = NULL;
   m_upperBestFitness = NULL;
   m_lowerBestFitness = NULL;
   m_tmpDivisionSize = NULL;
   m_increment = NULL;
} /* end PAES_AdaptiveGrid::AdaptiveGrid() */

/******************************************************************************
PAES_AdaptiveGrid::CTOR
******************************************************************************/
PAES_AdaptiveGrid::PAES_AdaptiveGrid(int depth, int numberOfFunctions) 
{                   
   m_numberOfFunctions  = numberOfFunctions;
   m_depth              = depth;
  
   m_currentGridSize = (long) floor(pow(2.0, (double)m_depth * (double)m_numberOfFunctions));
   m_hypercube    = new int[m_currentGridSize];
   m_divisionSize = new double[m_numberOfFunctions];
   m_gridLimits   = new double[m_numberOfFunctions];
  
   m_upperBestFitness = new double[m_numberOfFunctions];
   m_lowerBestFitness = new double[m_numberOfFunctions];  
  
   m_tmpDivisionSize = new double[m_numberOfFunctions];
   m_increment       = new int[m_numberOfFunctions];
  
   if (!m_currentGridSize  || !m_hypercube        || !m_divisionSize    || !m_gridLimits  || 
       !m_upperBestFitness || !m_lowerBestFitness || !m_tmpDivisionSize || !m_increment) 
   {
      printf("AdaptiveGrid::AdaptiveGrid-> Error when asking for memory\n");
      exit(-1);
   } /* end if() */
} /* end PAES_AdaptiveGrid::AdaptiveGrid() */

/******************************************************************************
PAES_AdaptiveGrid::DTOR
******************************************************************************/
PAES_AdaptiveGrid::~PAES_AdaptiveGrid(void) 
{
   delete [] m_hypercube;
   delete [] m_divisionSize;
   delete [] m_gridLimits;
   delete [] m_upperBestFitness;
   delete [] m_lowerBestFitness;
   delete [] m_tmpDivisionSize;
   delete [] m_increment;
} /* end PAES_AdaptiveGrid::~AdaptiveGrid() */

/******************************************************************************
PAES_AdaptiveGrid::updateGridLocations()
******************************************************************************/
void PAES_AdaptiveGrid::updateGridLocations(PAES_Population * population, PAES_Individual * individual) 
{
   int i;
   int j;
  
   for (i = 0; i < m_numberOfFunctions; ++i) 
   {
      m_upperBestFitness[i] = PAES_MIN_INT;
      m_lowerBestFitness[i] = PAES_MAX_INT;
   } /* end for() */
  
   for (i = 0; i < m_numberOfFunctions; i++) 
   {
      if (individual->getFitness()[i] < m_lowerBestFitness[i])
      {
         m_lowerBestFitness[i] = individual->getFitness()[i];
      }
      if (individual->getFitness()[i] > m_upperBestFitness[i])
      {
         m_upperBestFitness[i] = individual->getFitness()[i];
      }

      for (j = 0; j < population->getPopulationSize(); j ++) 
      {
         if (population->getIth(j)->getFitness()[i] < m_lowerBestFitness[i]) 
         {
            m_lowerBestFitness[i] = population->getIth(j)->getFitness()[i];
         }
         if (population->getIth(j)->getFitness()[i] > m_upperBestFitness[i]) 
         {
            m_upperBestFitness[i] = population->getIth(j)->getFitness()[i];      
         }
      }/* end for() */
    
      m_divisionSize[i] = (m_upperBestFitness[i] - m_lowerBestFitness[i]);
   } /* end for() */
  
   for (i = 0; i < (int)pow(2.0, m_numberOfFunctions * m_depth); i++)
   {
      m_hypercube[i] = 0;
   }
    
   int location;  
   m_mostCrowdedHypercube = 0;
   location = findLocation(individual);
   individual->m_gridLocation = location;
   m_hypercube[location]++;
   for (i = 0; i < population->getPopulationSize(); i++) 
   {
      location = findLocation(population->getIth(i));
      (population->getIth(i))->m_gridLocation = location;
      m_hypercube[location]++;
      if (m_hypercube[location] > m_hypercube[m_mostCrowdedHypercube])
      {
         m_mostCrowdedHypercube = location;
      }
   } /* end for() */
} /* end PAES_AdaptiveGrid::updateGridLocations() */

/******************************************************************************
PAES_AdaptiveGrid::findLocation()

Find the location of the individual in the adaptive grid
******************************************************************************/
int PAES_AdaptiveGrid::findLocation(PAES_Individual * individual) 
{
   int i;
   int j;
   int location;
   int counter;
  
   location = 0;
   counter  = 1;
   for (i = 0; i < m_numberOfFunctions; i++) 
   {
      m_increment[i] = counter;
      counter *= 2;
      m_tmpDivisionSize[i] = m_divisionSize[i];
   } /* end for() */
  
   for (i = 1; i <= m_depth; i++) 
   {
      for (j = 0; j < m_numberOfFunctions; j++) 
      {
         if (individual->getFitness()[j] < (m_tmpDivisionSize[j]/2 + m_lowerBestFitness[j]))
         {
            location = m_increment[j];
         }
         else
         {
            m_lowerBestFitness[j] += m_tmpDivisionSize[j] / 2.0;
         }
      }/* end for() */
      for (j = 0; j < m_numberOfFunctions; j++) 
      {
         m_increment[j] *= m_numberOfFunctions * 2;
         m_tmpDivisionSize[j] /= 2.0;
      } /* for() */
   } /* end for() */
  
   return location;
} /* end PAES_AdaptiveGrid::findLocation() */
        
