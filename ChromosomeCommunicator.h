/******************************************************************************
File     : ChromosomeCommunicator.h
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
08-17-04    lsm   RAM fragmentation fixes
03-21-05    lsm   Added support for user-defined seeding of initial population
01-01-07    lsm   Algorithms now use an abstract model base class (ModelABC).
******************************************************************************/
#ifndef CHROMOSOME_COMMUNICATOR_H
#define CHROMOSOME_COMMUNICATOR_H

#include "MyHeaderInc.h"

//forward declarations
class ModelABC;
class Chromosome;
class ParameterGroup;
class ModelBackup;

/******************************************************************************
class ChromosomeCommunicator

  Abstract class that acts as a link to a chromosome object.
******************************************************************************/
class ChromosomeCommunicator
{
   public:
      virtual ~ChromosomeCommunicator(void) { DBG_PRINT("ChromosomeCommunicator::DTOR"); }
      virtual void Destroy(void) = 0;
      virtual void EvalFitness(Chromosome * pChromo)=0;
      virtual Chromosome * CreateProto(double rate)=0;
      virtual ParameterGroup * ConvertChromosome(Chromosome * pChromo) = 0;
      virtual ParameterGroup * GetParamGroupPtr(void) = 0;
      virtual void MakeParameterCorrections(Chromosome * pChromo) = 0;
}; /* end class ChromosomeCommunicator */

/******************************************************************************
class ModelChromosomeCommunicator

  Class that has knowledge of both the model class and the chromosome class. 
  The class also understands the working of these details.
  
  For example - This class has the knowledge of the parameters and their type
  This class can be modified to be less dependent on the details of the 
  parameters etc but the for the sake of ease of implementation, this 
  dependency has been used. 
******************************************************************************/
class ModelChromoComm : public ChromosomeCommunicator
{
   public:
      void Destroy(void);
      ~ModelChromoComm(void){ DBG_PRINT("ModelChromoComm::DTOR"); Destroy(); }
      ModelChromoComm(ModelABC * pModel);
      void EvalFitness(Chromosome * pChromo);
      Chromosome * CreateProto(double rate);
      ParameterGroup * ConvertChromosome(Chromosome * pChromo);
	  ParameterGroup * GetParamGroupPtr(void);
      void SetMaxEvals(int maxEvals){ m_MaxEvals = maxEvals;}
      void MakeParameterCorrections(Chromosome * pChromo);

   private :
      ModelABC * m_pModel; 
      double * m_xb; //best parameter set
      int m_MaxEvals;        
}; /* end class ModelChromoComm */

#endif /* CHROMOSOME_COMMUNICATOR_H */

