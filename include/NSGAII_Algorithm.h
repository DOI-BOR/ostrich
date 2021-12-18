/******************************************************************************
File     : NSGAII_Algorithm.h
Author   : L. Shawn Matott and Kalyanmoy Deb
Copyright: 2018, L. Shawn Matott (adapted into OSTRICH)
           ????, Kalyanmoy Deb (original standalone implementation)

Non-dominated Sorting Genetic Algorithm - version 2.

Version History
01-24-18    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef NSGAII_H
#define NSGAII_H

#define NSGAII_SQUARE(x) ((x)*(x))

// parent class
#include "AlgorithmABC.h"

// forward decs
class ModelABC;

/******************************************************************************
NSGAII_Individual

A structure defining an individual
******************************************************************************/
typedef struct       
{
   int * genes; /* [maxchrom] bianry chromosome */
   float * gener; /* [maxchrom] real chromosomes */
   float dumfit; /* dummy fitness */
   int rank; /* Rank of the individual*/
   int flag; /* Flag for ranking */
   float * xreal; /* [maxvar] list of real variables */
   float * xbin; /* [maxvar] list of decoded value of the chromosome */
   float * fitness; /* [maxfun] Fitness values */
   float * constr; /* [maxcons] Constraints values */
   float cub_len; /* crowding distance of the individual */
   float error; /* overall constraint violation for the individual */
}NSGAII_Individual; 

/******************************************************************************
NSGAII_Population

A structure defining a population of individuals
******************************************************************************/
typedef struct
{
   int maxrank; /* Maximum rank present in the population */
   float * rankrat; /* [maxpop] Rank Ratio */
   int * rankno; /* [maxpop] Individual at different ranks */
   NSGAII_Individual * ind; /* [maxpop] Different Individuals */
   NSGAII_Individual * ind_ptr; 
}NSGAII_Population;

/******************************************************************************
NSGAII_GlobPop

Population structure for the pool having both the old as well as new population
******************************************************************************/
typedef struct
{
  int maxrank; /* Max rank of the global population */
  int ** rankar; /* [2*maxpop][2*maxpop] record of array of individual numbers at a particular rank */
  int * rankno; /* [2*maxpop] record of no. of individuals at a particular rank */
  int ** genes; /* [2*maxpop][maxchrom] */      
  int * rank; /* [2*maxpop] rank of different individuals */
  int * flag; /* [2*maxpop] Setting the flag */

  float ** fitness; /* [2*maxpop][maxfun] Fitness function values for the different individuals */
  float * cub_len; /* [2*maxpop] Dummyfitness */
  float ** xreal; /* [2*maxpop][maxvar] value of the decoded variables for different individuals */
  float ** xbin; /* [2*maxpop][maxvar] binray-coded variables */
  float * error; /* [2*maxpop] Error Values of the individuals */
  float ** constr; /* [2*maxpop][maxcons] */
}NSGAII_GlobPop;

/******************************************************************************
class NSGAII
******************************************************************************/
class NSGAII : public AlgorithmABC
{
   public:
      NSGAII(ModelABC * pModel);
      ~NSGAII(void){ DBG_PRINT("NSGAII::DTOR"); Destroy(); }
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);
      void WarmStart(void){ return;}
      int  GetCurrentIteration(void) { return m_CurIter; }

   private:
      void crossover(NSGAII_Population *new_pop_ptr, NSGAII_Population *mate_pop_ptr);
      void decode(NSGAII_Population *pop_ptr);
      void dumfitness(NSGAII_Population *pop_ptr);
      void share1(NSGAII_Population *pop_ptr,int rnk);
      void sort1(int m1);
      void func(NSGAII_Population *pop_ptr);
      void func_parallel(NSGAII_Population *pop_ptr, int myrank, int nprocs);
      void init(NSGAII_Population *pop_ptr);
      void keepalive(NSGAII_Population *pop1_ptr, NSGAII_Population *pop2_ptr, NSGAII_Population *pop3_ptr, int gen);
      void grank(int gen);
      void grankc(int gen);
      int indcmp1(float *ptr1,float *ptr2);
      void gsort(int rnk,int sel);
      void gshare(int rnk);
      void sort(int m1);
      void mutate(NSGAII_Population *new_pop_ptr);             
      void rankcon(NSGAII_Population *pop_ptr);
      int indcmp3(float *ptr1,float *ptr2);
      int flip(float prob);
      double noise(double mu , double sigma);
      double randomnormaldeviate(void);
      float randomperc(void); 
      int rnd(int low, int high); 
      float rndreal(float lo , float hi); 
      void ranking(NSGAII_Population *pop_ptr);
      int indcmp(float *ptr1,float *ptr2);
      void realcross(NSGAII_Population *new_pop_ptr, NSGAII_Population *mate_pop_ptr);
      void realinit(NSGAII_Population *pop_ptr);
      void real_mutate(NSGAII_Population *new_pop_ptr);             
      void rselect(NSGAII_Population *old_pop_ptr, NSGAII_Population *mate_pop_ptr);
      void report(int t, NSGAII_Population *pop1_ptr, NSGAII_Population *pop2_ptr,FILE *rep_ptr,FILE *gen_ptr, FILE *lastit);
      void roulette(NSGAII_Population *old_pop_ptr, NSGAII_Population *mate_pop_ptr);
      void nselect(NSGAII_Population *old_pop_ptr, NSGAII_Population *pop2_ptr);
      void unicross(NSGAII_Population *new_pop_ptr, NSGAII_Population *mate_pop_ptr);

      void Delete_NSGAII_GlobPop(NSGAII_GlobPop * pGlobalPop);
      void Delete_NSGAII_Population(NSGAII_Population * pPop);
      void Delete_NSGAII_Individual(NSGAII_Individual * pInd);

      void Create_NSGAII_GlobPop(NSGAII_GlobPop * pGlobalPop);
      void Create_NSGAII_Population(NSGAII_Population * pPop);
      void Create_NSGAII_Individual(NSGAII_Individual * pInd);

      int UpdateLists(double * pX, int nX, double * pF, int nF);

      void EvaluateSamples(int nSamples);

      int m_maxpop;  /* Max population */
      int m_maxchrom;  /* Max chromosome length (bits) */
      int m_maxvar;  /* Max no. of variables */
      int m_maxfun;  /* Max no. of functions */
      int m_maxcons;  /* Max no. of Constraints */
      int m_gener; /* No of generations */
      int m_nvar; /* No of variables */
      int m_nchrom; /* No of chromosomes */
      int m_ncons; /* No of Constraints */
      int * m_vlen; /* [maxvar]  Array to store no of bits for each variable */
      int m_nmut; /* No of Mutations */
      int m_ncross; /* No of crossovers */
      int m_ans; 

      float m_seed; /* Random Seed */
      float m_pcross; /* Cross-over Probability */
      float m_pmut_b; /* Mutation Probability - binary vars */
      float m_pmut_r; /* Mutation Probability - real vars */
      float ** m_lim_b; /* [maxvar][2] limits of binary vars */
      float ** m_lim_r; /* [maxvar][2] limits of real vars */
      float m_di; /* Distribution Index for the Cross-over */
      float m_dim; /* Distribution Index for the Mutation */
      float m_delta_fit; /* variables required for fitness sharing */
      float m_min_fit;
      float m_front_ratio;

      int m_optype; /* Cross-over type */
      int m_nfunc;  /* No of functions */
      int m_sharespace; /* Sharing space (either parameter or fitness) */

      double * m_coef; /* [maxvar] Variable used for decoding */

      int m_popsize; /* Population Size */
      int m_chrom; /*Chromosome size (bits)*/

      /* various populations */
      NSGAII_Population m_oldpop;
      NSGAII_Population m_newpop;
      NSGAII_Population m_matepop;
      NSGAII_Population * m_old_pop_ptr;
      NSGAII_Population * m_new_pop_ptr;
      NSGAII_Population * m_mate_pop_ptr;

      int ** m_rec_arr;  /* [maxpop][maxpop] This array holds the record of indiviuals at different ranks*/
      int * m_r_arr; /* [maxpop] This array gives the number of individuals at different ranks*/
      float ** m_fpara; /* [maxpop][2] Stores individual no and its fitness in one dimension*/

      NSGAII_GlobPop m_globalpop;
      NSGAII_GlobPop *m_global_pop_ptr;

      int m_left;
      int m_Lastrank;
      float ** m_fpara1; /* [2*maxpop][2] */

      FILE * m_rep_ptr;

      ModelABC * m_pModel;
      int m_DomListSize; //number of dominated solutions
      int m_NonDomListSize; //number of non-dominated solutions
      ArchiveStruct * m_pNonDomList; //non-dominated solutions
      ArchiveStruct * m_pDomList; //dominated solutions
      int m_CurIter;
}; /* end class NSGAII */

extern "C" {
void NSGAII_Program(int argC, StringType argV[]);
}

#endif /* NSGAII_H */
