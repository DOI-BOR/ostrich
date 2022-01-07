/******************************************************************************
File     : NSGAII_Algorithm.h
Author   : L. Shawn Matott and Kalyanmoy Deb
Copyright: 2018, L. Shawn Matott (adapted into OSTRICH)
           ????, Kalyanmoy Deb (original standalone implementation)

Non-dominated Sorting Genetic Algorithm - version 2.

The program generates a few output files. These are described as
1.output.out
   This file has the detailed record for all the variables,
   the fitness values, constraint values, overall constraint
   violation (penalty)  and their ranks for all the members
   of old population in the left hand side of the |**|
   and of new population in the right hand side.

2.all_fitness.out
   This file prints the record of all the fitness values for
   different individual of new popultion created at all
   generations.

3.g_rank_record.out
   This file maintains the record of individuals in global pop-
   -ulation at different ranks for all the generations.

4.ranks.out
   This file prints the number of individual at different ranks
   in old and new population and finds rank ratios

5.final_fitness.out
   This file has the fitness value of all feasible and
   non-dominated individuals at the final generation

6.final_var.out
   This file has the all the variables of the feasible
   and non-dominated individuals at the final generation.
   The i-th solutions here corresponds to the i-th solution
   in the final_fitness.out file. 

7.plot.out       
   This file contains gnuplot-based file for plotting
   the non-dominated feasible solutions obtained by the code.

   It is recommended to delete or rename all the *.out files
   obtained from the previous runs as some files are opened in
   append mode so they give false resemblence of data if the
   user is not careful

Version History
01-24-18    lsm   added copyright information and initial comments.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "NSGAII_Algorithm.h"
#include "Model.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
NSGAII::NSGAII(ModelABC * pModel)
{
   int i;
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pNonDomList = NULL;
   m_pDomList = NULL;
   m_NonDomListSize = 0;
   m_DomListSize = 0;
   m_CurIter = 0;

   m_maxpop = 500;  /* Max population */
   m_maxchrom = 0;  /* Max chromosome length (bits) */
   m_maxvar = pModel->GetParamGroupPtr()->GetNumParams();  /* Max no. of variables */
   m_maxfun = pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);  /* Max no. of functions */
   m_maxcons = 0;  /* Max no. of Constraints */
   m_gener = 50; /* No of generations */
   m_nvar = m_maxvar; /* No of variables (i.e. real coded parameters) */
   m_nchrom = 0; /* No of chromosomes (i.e. binary coded parameters) */
   m_ncons = m_maxcons; /* No of Constraints */
   m_vlen = new int[m_maxvar]; /* [maxvar]  Array to store no of bits for each variable */
   m_nmut = 0; /* No of Mutations */
   m_ncross = 0; /* No of crossovers */
   m_ans = 1; /* alwys use limits */ 
   m_seed = 0; /* Random Seed */
   m_pcross = (float)0.1; /* Cross-over Probability */
   m_pmut_b = (float)0.05; /* Mutation Probability - binary vars */
   m_pmut_r = (float)0.05; /* Mutation Probability - real vars */
   /* assign parameter bounds */
   m_lim_b = new float *[m_maxvar];
   m_lim_r = new float *[m_maxvar];
   for(i = 0; i < m_maxvar; i++)
   {
      m_lim_b[i] = new float[2];
      m_lim_r[i] = new float[2];
      m_lim_b[i][0] = (float)(pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLwrBnd());
      m_lim_b[i][1] = (float)(pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUprBnd());
      m_lim_r[i][0] = m_lim_b[i][0];
      m_lim_r[i][1] = m_lim_b[i][1];
      m_vlen[i] = 30; /* default to 30 bits per chromosome */
   }
   m_di = 20.00; /* Distribution Index for the Cross-over */
   m_dim = 10.00; /* Distribution Index for the Mutation */
   m_delta_fit = 0.00; /* variables required for fitness sharing */
   m_min_fit = 0.00;
   m_front_ratio = 0.00;

   m_optype = 1; /* Cross-over type */
   m_nfunc = m_maxfun;  /* No of functions */
   m_sharespace = 0; /* Sharing space (either parameter or fitness) */

   m_coef = new double[m_maxvar]; /* [maxvar] Variable used for decoding */

   m_popsize = m_maxpop; /* Population Size */
   m_chrom = 0; /*Chromosome size (bits) */

   /* various populations */
   m_old_pop_ptr = NULL;
   m_new_pop_ptr = NULL;
   m_mate_pop_ptr = NULL;

   /* defer allocations untile input file is processed ... */
   m_rec_arr = NULL;  /* [maxpop][maxpop] This array holds the record of indiviuals at different ranks*/
   m_r_arr = NULL; /* [maxpop] This array gives the number of individuals at different ranks*/
   m_fpara = NULL; /* [maxpop][2] Stores individual no and its fitness in one dimension*/
   m_fpara1 = NULL; /* [2*maxpop][2] */

   m_global_pop_ptr = NULL;
   m_rep_ptr = NULL;

   m_left = 0;
   m_Lastrank = 0;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by NSGAII and it's member variables.
******************************************************************************/
void NSGAII::Destroy(void)
{
   int i;
   ArchiveStruct * pCur, * pDel;

   for(pCur = m_pNonDomList; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   for(pCur = m_pDomList; pCur != NULL;)
   {
     pDel = pCur;
     pCur = pCur->pNext;
     delete [] pDel->F;
     delete [] pDel->X;
     delete pDel; 
   }

   delete [] m_vlen; /* [maxvar]  Array to store no of bits for each variable */

   /* parameter bounds */
   for(i = 0; i < m_maxvar; i++)
   {
      delete [] m_lim_b[i];
      delete [] m_lim_r[i];
   }
   delete [] m_lim_b;
   delete [] m_lim_r;

   delete [] m_coef; /* [maxvar] Variable used for decoding */

   /* various populations: */
   Delete_NSGAII_Population(&m_oldpop);
   Delete_NSGAII_Population(&m_newpop);
   Delete_NSGAII_Population(&m_matepop);

   for(i = 0; i < m_maxpop; i++)
   {
      delete [] m_rec_arr[i];
      delete [] m_fpara[i];
   }
   delete [] m_rec_arr; /* [maxpop][maxpop] This array holds the record of indiviuals at different ranks*/
   delete [] m_r_arr; /* [maxpop] This array gives the number of individuals at different ranks*/
   delete [] m_fpara;; /* [maxpop][2] Stores individual no and its fitness in one dimension*/

   for(i = 0; i < 2*m_maxpop; i++)
   {
      delete [] m_fpara1[i]; /* [2*maxpop][2] */
   }
   delete [] m_fpara1;
 
   /* use helper function for cleaning up this data struct */
   Delete_NSGAII_GlobPop(&m_globalpop);;

   if(m_rep_ptr != NULL)
   {
      fclose(m_rep_ptr);
   }

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Delete_NSGAII_Individual()

Free up memory in population struct.
******************************************************************************/
void NSGAII::Delete_NSGAII_Individual(NSGAII_Individual * pInd)
{
   delete [] pInd->genes; /* [maxchrom] bianry chromosome */
   delete [] pInd->gener; /* [maxchrom] real chromosomes */
   delete [] pInd->xreal; /* [maxvar] list of real variables */
   delete [] pInd->xbin; /* [maxvar] list of decoded value of the chromosome */
   delete [] pInd->fitness; /* [maxfun] Fitness values */
   delete [] pInd->constr; /* [maxcons] Constraints values */
}/* end Delete_NSGAII_Individual() */

/******************************************************************************
Delete_NSGAII_Population()

Free up memory in population struct.
******************************************************************************/
void NSGAII::Delete_NSGAII_Population(NSGAII_Population * pPop)
{
   int i;
   delete [] pPop->rankrat; /* [maxpop] Rank Ratio */
   delete [] pPop->rankno; /* [maxpop] Individual at different ranks */
   for(i = 0; i < m_maxpop; i++)
   {
      Delete_NSGAII_Individual(&(pPop->ind[i])); /* [maxpop] Different Individuals */
   }
}/* end Delete_NSGAII_Population() */


/******************************************************************************
UpdateLists()

Update the dominated and non-dominated lists with latest sample.
******************************************************************************/
int NSGAII::UpdateLists(double * pX, int nX, double * pF, int nF)
{
   int i;
   double Fcur, Ftst;
   ArchiveStruct * pArch, *pCur, *pPrev, *pNxt;
   bool bDominates, bIsDominated, bMarkForInsertion;
   pArch = new ArchiveStruct;
   pArch->F = pF;
   pArch->X = pX;
   pArch->nX = nX;
   pArch->nF = nF;
   pArch->pNext = NULL;

   //first entry is always non-dominated
   if ((m_DomListSize == 0) && (m_NonDomListSize == 0))
   {
      m_pDomList = NULL;
      m_pNonDomList = pArch;
      m_NonDomListSize++;
      return ARCHIVE_NON_DOM;
   }

   //assume solution is non-dominated until we discover otherwise
   bMarkForInsertion = true;

   //compare against current list of non-dominated solutions
   ArchiveStruct * pDummy;
   for (pCur = m_pNonDomList; pCur != NULL;)
   {
      //save next item since pCur->pNext may be changed during processing
      pDummy = pCur->pNext;

      //does new solution (Ftst) dominate the existing solution (Fcur)?
      bDominates = true;
      for (i = 0; i < pArch->nF; i++)
      {
         Fcur = pCur->F[i];
         Ftst = pArch->F[i];
         if (Fcur < Ftst)
         {
            bDominates = false;
            break;
         }/* end if() */
      }/* end for() */

      //is new solution (Ftst) dominated by an existing solution (Fcur)?
      if (bDominates == false)
      {
         bIsDominated = true;
         for (i = 0; i < pArch->nF; i++)
         {
            Fcur = pCur->F[i];
            Ftst = pArch->F[i];
            if (Ftst < Fcur)
            {
               bIsDominated = false;
               break;
            }/* end if() */
         }/* end for() */
      }/* end if() */

      /* -----------------------------------------------------------------------------
      Existing solution is dominated. Remove it from list of non-dominated solutions
      and mark new solution for insertion into the non-dominated list.
      ----------------------------------------------------------------------------- */
      if (bDominates == true)
      {
         //solution to be removed is at head of list.
         if (pCur == m_pNonDomList)
         {
            m_pNonDomList = pCur->pNext;
            m_NonDomListSize--;
         }/* end if() */
         else //somewhere in middle or at end.
         {
            for (pPrev = m_pNonDomList; pPrev->pNext != pCur;)
            {
               pPrev = pPrev->pNext;
            }
            pNxt = pCur->pNext;
            pPrev->pNext = pNxt; //removes item from list
            m_NonDomListSize--;
         }/* end else() */

         //insert at head of dominated list
         pCur->pNext = NULL;
         pNxt = m_pDomList;
         m_pDomList = pCur;
         m_pDomList->pNext = pNxt;
         m_DomListSize++;
      }/* end if() */
      /* -----------------------------------------------------------------------------
      New solution is dominated. Make note so that it is not inserted into the list
      of non-dominated solutions.
      ----------------------------------------------------------------------------- */
      else if (bIsDominated == true)
      {
         bMarkForInsertion = false;
      }/* end if() */

      //advance to next item
      pCur = pDummy;
   }/* end for() */

   //insert new solution into list of non-dominated solutions?
   if (bMarkForInsertion == true)
   {
      //insert at head of non-dominated list
      pNxt = m_pNonDomList;
      m_pNonDomList = pArch;
      m_pNonDomList->pNext = pNxt;
      m_NonDomListSize++;
      return ARCHIVE_NON_DOM;
   }/* end if() */
   else
   {
      //insert at head of dominated list
      pNxt = m_pDomList;
      m_pDomList = pArch;
      m_pDomList->pNext = pNxt;
      m_DomListSize++;
      return ARCHIVE_DOM;
   }/* end else() */
}/* end UpdateLists() */

/******************************************************************************
Calibrate()

Solve the multi-criteria Least-Squares minimization problem using NSGAII.
******************************************************************************/
void NSGAII::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Search for the pareto front using NSGAII.
******************************************************************************/
void NSGAII::Optimize(void)
{
   int num, nobj;
   StatusStruct pStatus;
   ParameterGroup * pGroup;

   /*Some Local variables to this Problem (Counters And some other pointers*/
   int i,j,l,f,maxrank1;
   //float * ptr;
   float tot;
   /*File Pointers*/
   FILE * gen_ptr;
   FILE * rep2_ptr;
   FILE * end_ptr;
   FILE * g_var;
   FILE * lastit;

   InitFromFile(GetInFileName());
   
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();
   nobj = m_pModel->GetObjFuncPtr()->CalcMultiObjFunc(NULL, -1);

   WriteSetup(m_pModel, "NSGAII - Non-dominated Sorted Genatic Algorithm - version 2");
   //write banner
   WriteBanner(m_pModel, "gen   ", "Convergence Value");

   /*Opening various NSGAII output files*/
   gen_ptr =fopen("NSGAII_all_fitness.out","w");
   rep2_ptr = fopen("NSGAII_ranks.out","w");
   end_ptr = fopen("NSGAII_final_fitness.out","w");
   g_var = fopen("NSGAII_final_var.out","w");
   lastit = fopen("NSGAII_plot.out","w");

   m_old_pop_ptr = &(m_oldpop);

   m_nmut = 0;
   m_ncross = 0;

   fprintf(m_rep_ptr,"Results in a file\n");
   fprintf(end_ptr,"# Last generation population (Feasible and non-dominated)\n");
   fprintf(end_ptr,"# Fitness_vector (first %d)  Constraint_violation (next %d)  Overall_penalty\n",m_nfunc,m_ncons);
   fprintf(g_var,"#Feasible Variable_vectors for non-dominated solutions at last generation\n");
   fprintf(g_var,"# Real (first %d)  Binary (next %d)\n",m_nvar,m_nchrom);
   fprintf(lastit,"# Feasible and Non-dominated Objective Vector\n");

   /*Binary Initializaton*/
   if(m_nchrom > 0)
   {
      init(m_old_pop_ptr);  
   }
   if(m_nvar > 0)
   {
      realinit(m_old_pop_ptr);
   } 
   m_old_pop_ptr = &(m_oldpop);

   // decode binary strings
   decode(m_old_pop_ptr); 

   m_old_pop_ptr = &(m_oldpop);
   m_new_pop_ptr = &(m_newpop);
  
   for(j = 0; j < m_popsize; j++)
   {
      /*-------------------------------------------------------
      Initializing the Rank array having different individuals
      at a particular rank to zero
      -------------------------------------------------------*/
      m_old_pop_ptr->rankno[j] = 0;
      m_new_pop_ptr->rankno[j] = 0;
   }
  
   m_old_pop_ptr = &(m_oldpop);

   /* Function Calculation */  
   WriteInnerEval(WRITE_GA, m_popsize, '.');

   func(m_old_pop_ptr); /* func() evaluates entire population */

   /* update OSTRICH run record and status file */
   i = 0;
   WriteInnerEval(WRITE_ENDED, 0, '.');
   pStatus.pct = ((float)100.00*(float)(i)) / (float)m_gener;
   pStatus.numRuns = (i + 1)*m_popsize;
   WriteMultiObjRecord(m_pModel, i, m_pNonDomList, pStatus.pct);
   WriteStatus(&pStatus);

   fprintf(m_rep_ptr,"----------------------------------------------------\n");
   fprintf(m_rep_ptr,"Statistics at Generation 0 ->\n");
   fprintf(m_rep_ptr,"--------------------------------------------------\n");

   //main optimization loop   
   /********************************************************************/
   /*----------------------GENERATION STARTS HERE----------------------*/
   for (i = 0; i < m_gener; i++)
   {
      pStatus.curIter = m_CurIter = i+1;
      if(IsQuit() == true){ break;}

      m_old_pop_ptr = &(m_oldpop);
      m_mate_pop_ptr = &(m_matepop);
      fprintf(m_rep_ptr,"Population at generation no. -->%d\n",i+1);
      fprintf(gen_ptr,"#Generation No. -->%d\n",i+1);
      fprintf(gen_ptr,"#Variable_vector  Fitness_vector Constraint_violation Overall_penalty\n");
      
      /*--------SELECT----------------*/
      nselect(m_old_pop_ptr, m_mate_pop_ptr );
      
      m_new_pop_ptr = &(m_newpop);
      m_mate_pop_ptr = &(m_matepop);
      
      /*CROSSOVER----------------------------*/      
      if(m_nchrom > 0) 
      {
         if(m_optype == 1)
         {
            crossover(m_new_pop_ptr, m_mate_pop_ptr); /*Binary Cross-over*/
         }
	  
         if(m_optype == 2)
         {
            unicross(m_new_pop_ptr, m_mate_pop_ptr); /*Binary Uniform Cross-over*/
         }
      }/* end if() */
      if(m_nvar > 0) 
      {
         realcross(m_new_pop_ptr, m_mate_pop_ptr); /*Real Cross-over*/
      }
      
      /*------MUTATION-------------------*/
      m_new_pop_ptr = &(m_newpop);
      
      if(m_nchrom > 0)
      {
         mutate(m_new_pop_ptr); /*Binary Mutation */
      }

      if(m_nvar > 0)
      {
         real_mutate(m_new_pop_ptr); /*Real Mutation*/
      }

      m_new_pop_ptr = &(m_newpop);
      
      /*-------DECODING----------*/
      if(m_nchrom > 0)
      {
         decode(m_new_pop_ptr); /*Decoding for binary strings*/
      }
      
      /*----------FUNCTION EVALUATION-----------*/
      WriteInnerEval(WRITE_GA, m_popsize, '.');

      m_new_pop_ptr = &(m_newpop);
      func(m_new_pop_ptr);
      
      /*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
      m_old_pop_ptr = &(m_oldpop);
      m_new_pop_ptr = &(m_newpop);
      m_mate_pop_ptr = &(m_matepop);
      
      /*Elitism And Sharing Implemented*/
      keepalive(m_old_pop_ptr,m_new_pop_ptr,m_mate_pop_ptr,i+1);      
      
      m_mate_pop_ptr = &(m_matepop);
      if(m_nchrom > 0)
      {
         decode(m_mate_pop_ptr);
      }
      m_mate_pop_ptr = &(m_matepop);

      /*------------------REPORT PRINTING--------------------------------*/  
      report(i, m_old_pop_ptr, m_mate_pop_ptr, m_rep_ptr, gen_ptr, lastit);
      
      /*==================================================================*/
      
      /*----------------Rank Ratio Calculation------------------------*/
      m_new_pop_ptr = &(m_matepop);
      m_old_pop_ptr = &(m_oldpop);
      
      /*Finding the greater maxrank among the two populations*/
      if(m_old_pop_ptr->maxrank > m_new_pop_ptr->maxrank)
      {
         maxrank1 = m_old_pop_ptr->maxrank;
      }
      else 
      {
         maxrank1 = m_new_pop_ptr->maxrank;
      }
      fprintf(rep2_ptr,"--------RANK AT GENERATION %d--------------\n",i+1);
      fprintf(rep2_ptr,"Rank old ranks   new ranks     rankratio\n");

      for(j = 0; j < maxrank1; j++)
      { 
         /*Sum of the no of individuals at any rank in old population and the new populaion*/
         tot = (float)((m_old_pop_ptr->rankno[j])+ (m_new_pop_ptr->rankno[j]));
	  
         /*Finding the rank ratio for new population at this rank*/
         m_new_pop_ptr->rankrat[j] = (m_new_pop_ptr->rankno[j])/tot;
	  
         /*Printing this rank ratio to a file called ranks.dat*/
         fprintf(rep2_ptr," %d\t  %d\t\t %d\t %f\n",j+1,m_old_pop_ptr->rankno[j],m_new_pop_ptr->rankno[j],m_new_pop_ptr->rankrat[j]);	  
      }
      
      fprintf(rep2_ptr,"-----------------Rank Ratio-------------------\n");
      /*==================================================================*/
      
      /*=======Copying the new population to old population======*/
      m_old_pop_ptr = &(m_oldpop);
      m_new_pop_ptr = &(m_matepop);

      for(j = 0; j < m_popsize; j++)
      {
         m_old_pop_ptr->ind_ptr = &(m_old_pop_ptr->ind[j]);
         m_new_pop_ptr->ind_ptr = &(m_new_pop_ptr->ind[j]);
         if(m_nchrom > 0)
         {
            /*For Binary GA copying of the chromosome*/      
            for(l = 0; l < m_chrom; l++)
            {
               m_old_pop_ptr->ind_ptr->genes[l] = m_new_pop_ptr->ind_ptr->genes[l];
	         }     
            for(l = 0; l < m_nchrom; l++)
            {
               m_old_pop_ptr->ind_ptr->xbin[l] = m_new_pop_ptr->ind_ptr->xbin[l];
            }
         }/* end if() */
         if(m_nvar > 0)
         {
            /*For Real Coded GA copying of the chromosomes*/
            for(l = 0; l < m_nvar; l++)
            {
               m_old_pop_ptr->ind_ptr->xreal[l] = m_new_pop_ptr->ind_ptr->xreal[l];
            }
         }
	  
         /*Copying the fitness vector */	  
         for(l = 0 ; l < m_nfunc ;l++)
         {
            m_old_pop_ptr->ind_ptr->fitness[l] = m_new_pop_ptr->ind_ptr->fitness[l];
	      }
  
         /*Copying the dummy fitness*/
         m_old_pop_ptr->ind_ptr->cub_len = m_new_pop_ptr->ind_ptr->cub_len;
	  
         /*Copying the rank of the individuals*/
         m_old_pop_ptr->ind_ptr->rank = m_new_pop_ptr->ind_ptr->rank;
	  
         /*Copying the error and constraints of the individual*/
         m_old_pop_ptr->ind_ptr->error = m_new_pop_ptr->ind_ptr->error;

         for(l = 0;l < m_ncons;l++)
         {
            m_old_pop_ptr->ind_ptr->constr[l] = m_new_pop_ptr->ind_ptr->constr[l];
         }
	  
         /*Copying the flag of the individuals*/
         m_old_pop_ptr->ind_ptr->flag = m_new_pop_ptr->ind_ptr->flag;
      } /* end for(j) */
      
      maxrank1 = m_new_pop_ptr->maxrank;
	  
      /*Copying the array having the record of the individual at different ranks */
      for(l = 0;l < m_popsize;l++)
      {
         m_old_pop_ptr->rankno[l] = m_new_pop_ptr->rankno[l];
      }

      /*Copying the maxrank */
      m_old_pop_ptr->maxrank = m_new_pop_ptr->maxrank;
      
      /*Printing the fitness record for last generation in a file last*/
      if(i == (m_gener - 1)) // for the last generation
      {   
         m_old_pop_ptr = &(m_matepop);
         for(f = 0; f < m_popsize; f++) // for printing
         {
            m_old_pop_ptr->ind_ptr = &(m_old_pop_ptr->ind[f]);
	      
            // for all feasible solutions and non-dominated solutions
            if ((m_old_pop_ptr->ind_ptr->error <= 0.0) && (m_old_pop_ptr->ind_ptr->rank == 1))  
            {
               for(l = 0; l < m_nfunc; l++)
               {
                  fprintf(end_ptr,"%f\t",m_old_pop_ptr->ind_ptr->fitness[l]);
               }
               for(l = 0;l < m_ncons;l++)
               {
                  fprintf(end_ptr,"%f\t",m_old_pop_ptr->ind_ptr->constr[l]);
               }
               if (m_ncons > 0)
               {
                  fprintf(end_ptr,"%f\t",m_old_pop_ptr->ind_ptr->error);
               }
               fprintf(end_ptr,"\n");
		  
               if (m_nvar > 0)
               {
                  for(l = 0;l < m_nvar ;l++)
                  {
                     fprintf(g_var,"%f\t",m_old_pop_ptr->ind_ptr->xreal[l]);
                  }
                  fprintf(g_var,"  ");
               }
		  
               if(m_nchrom > 0)
               {
                  for(l = 0; l < m_nchrom; l++)
                  {
                     fprintf(g_var,"%f\t",m_old_pop_ptr->ind_ptr->xbin[l]);
                  }
               }/* end if() */
               fprintf(g_var,"\n");
            } /* end if (feasible) */
         } /* end for(f) */	  
      } /* end if(i == gener-1)) */

      /* update OSTRICH run record and status file */
      WriteInnerEval(WRITE_ENDED, (i+1), '.');

      pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_gener;
      pStatus.numRuns = (i+1)*m_popsize;
      WriteMultiObjRecord(m_pModel, (i+1), m_pNonDomList, pStatus.pct);
      WriteStatus(&pStatus);
   }/* end for(i) */
   /*                   Generation Loop Ends                                */
   /************************************************************************/
  
   WriteMultiObjOptimal(m_pModel, m_pNonDomList, m_pDomList);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   WriteAlgMetrics(this);
  
   fprintf(m_rep_ptr,"NO. OF CROSSOVER = %d\n",m_ncross);
   fprintf(m_rep_ptr,"NO. OF MUTATION = %d\n",m_nmut);
   fprintf(m_rep_ptr,"------------------------------------------------------------\n");
   fprintf(m_rep_ptr,"---------------------------------Thanks---------------------\n");
   fprintf(m_rep_ptr,"-------------------------------------------------------------\n");
   //printf("NOW YOU CAN LOOK IN THE FILE OUTPUT2.DAT\n");
  
   /*Closing the files*/
   fclose(m_rep_ptr);
   fclose(gen_ptr);
   fclose(rep2_ptr);
   fclose(end_ptr);
   fclose(g_var);
   fclose(lastit);
   m_rep_ptr = NULL;
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void NSGAII::WriteMetrics(FILE * pFile) 
{
   int i;

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : NSGAII - Non-dominated Sorted Genetic Algorithm - version 2\n");
   fprintf(pFile, "Num Generations         : %d\n", m_gener);
   fprintf(pFile, "Population Size         : %d\n", m_popsize);
   fprintf(pFile, "Non-Dominated Solutions : %d\n", m_NonDomListSize);  
   fprintf(pFile, "Dominated Solutions     : %d\n", m_DomListSize);     
   fprintf(pFile, "Crossover Dist. Index   : %f\n", m_di);
   fprintf(pFile, "Mutation Dist. Index    : %f\n", m_dim);
   fprintf(pFile, "Crossover Probability   : %f\n", m_pcross);
   fprintf(pFile, "Mutation Probability    : %f\n", m_pmut_r);
   fprintf(pFile, "Encoding Type           : ");
   
   if(m_nvar == m_maxvar)
   {
      fprintf(pFile, "real\n");
   }
   else
   {
      fprintf(pFile, "binary\n");
      fprintf(pFile, "Bits Per Parameter\n");
      for(i = 0; i < m_maxvar; i++)
      {
         fprintf(pFile, "   %s : %d\n", m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetName(), m_vlen[i]);
      }/* end for() */
      fprintf(pFile, "\n");
      fprintf(pFile, "Total Chromosome Length : %d\n", m_chrom);
   }
   fprintf(pFile, "Crossover Type          : ");
   if(m_optype == 1)
   {
      fprintf(pFile, "single-point\n");
   }
   else
   {
      fprintf(pFile, "uniform\n");
   }

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void NSGAII::InitFromFile(IroncladString pFileName)
{
   int i, j;
   FILE * pFile;
   char * line;
   char * pTmp;
   char tmp[DEF_STR_SZ];
   char encodingStr[DEF_STR_SZ];
   char crossoverStr[DEF_STR_SZ];
   
   //read in NSGAII configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open NSGAII config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginNSGAII", pFileName) == true)
   {
      FindToken(pFile, "EndNSGAII", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginNSGAII", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndNSGAII") == NULL)
      {
         if(strstr(line, "PopulationSize") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_popsize);
            m_maxpop = m_popsize;
         }
         else if(strstr(line, "NumGenerations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_gener);
         }
         else if(strstr(line, "CrossoverDistributionIndex") != NULL)
         {
            sscanf(line, "%f", &m_di);
         }
         else if(strstr(line, "MutationDistributionIndex") != NULL)
         {
            sscanf(line, "%f", &m_dim);
         }
         else if(strstr(line, "CrossoverProb") != NULL)
         {
            sscanf(line, "%s %f", tmp, &m_pcross);
         }
         else if(strstr(line, "MutationProb") != NULL)
         {
            sscanf(line, "%s %f", tmp, &m_pmut_r);
            m_pmut_b = m_pmut_r;
         }
         else if(strstr(line, "EncodingType") != NULL)
         {
            sscanf(line, "%s %s", tmp, encodingStr);
            MyStrLwr(encodingStr);
            if(strcmp(encodingStr, "real") == 0)
            {
               m_nvar = m_maxvar;
               m_nchrom = 0;
            } 
            else if(strcmp(encodingStr, "binary") == 0)
            {
               m_nvar = 0;
               m_nchrom = m_maxvar;
            }
            else /* default to real encoding */
            {
               m_nvar = m_maxvar;
               m_nchrom = 0;
            }
         }/*end else if() */         
         else if(strstr(line, "CrossoverType") != NULL)
         {
            sscanf(line, "%s %s", tmp, crossoverStr);
            MyStrLwr(crossoverStr);
            if(strcmp(crossoverStr, "single-point") == 0)
            {
               m_optype = 1;
            } 
            else if(strcmp(encodingStr, "uniform") == 0)
            {
               m_optype = 2;
            }
            else /* default to single-point */
            {
               m_optype = 1;
            }
         }/*end else if() */         
         else if(strstr(line, "BitsPerParameter") != NULL)
         {
            MyStrLwr(line);
            pTmp = line + (int)strlen("BitsPerParameter");

            for(i = 0; i < m_maxvar; i++)
            {
                j = ExtractString(pTmp, tmp);
                m_vlen[i] = atoi(tmp);
                if(j == -1){ break;}
                pTmp += j;               
            }

            for(i = (i+1); i < m_maxvar; i++)
            { 
               m_vlen[i] = m_vlen[0];
            }
         }/*end else if() */         
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   fclose(pFile);

   /*----------------------------------------------------------------
   For Binary Coded Problem --- Specify the no of bits for each 
   variable. Total Sum of the bit value must be equal to chromosome 
   length.
   -----------------------------------------------------------------*/
   m_chrom = 0;
   if (m_nchrom > 0)
   {      
      for (i = 0; i < m_nchrom; i++)
      {
         m_chrom += m_vlen[i];
	   }
   }
   m_maxchrom = m_chrom;

   /* Give the initial seed */
   m_seed = ((float)ReadRandomSeed() / (float)MY_RAND_MAX);

   /* allocate space for various populations */
   Create_NSGAII_Population(&m_oldpop);
   Create_NSGAII_Population(&m_newpop);
   Create_NSGAII_Population(&m_matepop);
   Create_NSGAII_GlobPop(&m_globalpop);
  
   /* other deferred allocations */
   m_r_arr = new int[m_maxpop]; /* [maxpop] This array gives the number of individuals at different ranks */
   m_rec_arr = new int *[m_maxpop]; /* [maxpop][maxpop] This array holds the record of indiviuals at different ranks */
   for(i = 0; i < m_maxpop; i++)
   {
      m_rec_arr[i] = new int[m_maxpop];
   }/* end for() */

   m_fpara = new float *[m_maxpop]; /* [maxpop][2] Stores individual no and its fitness in one dimension */
   for(i = 0; i < m_maxpop; i++)
   {
      m_fpara[i] = new float[2];
   }/* end for() */

   m_fpara1 = new float *[2*m_maxpop]; /* [2*maxpop][2] */
   for(i = 0; i < (2*m_maxpop); i++)
   {
      m_fpara1[i] = new float[2];
   }/* end for() */

   /*Print the GA parameters and problem parameters in the file output.dat*/
   m_rep_ptr = fopen("NSGAII_output.out","w");
   fprintf(m_rep_ptr,"GA PARAMETERS\n");
   fprintf(m_rep_ptr,"-------------------------------------------------------\n");
  
   fprintf(m_rep_ptr,"Population Size ->%d\n",m_popsize);
  
   fprintf(m_rep_ptr,"Chromosome Length ->%d\n",m_chrom);
  
   fprintf(m_rep_ptr,"No. of generations ->%d\n",m_gener);
 
   fprintf(m_rep_ptr,"No. of Functions ->%d\n",m_nfunc);

   fprintf(m_rep_ptr,"No. of Constraints ->%d\n",m_ncons);
  
   if(m_nchrom > 0)
   {
      fprintf(m_rep_ptr,"No. of binary-coded variables ->%d\n", m_nchrom);
   }
   if(m_nvar > 0) 
   {
      fprintf(m_rep_ptr,"No. of real-coded variables ->%d\n", m_nvar);
   }
   fprintf(m_rep_ptr,"Selection Strategy is Tournament Selection\n");
    
   for(i = 0; i < m_nchrom; i++)
   {
      fprintf(m_rep_ptr,"Binary-coded variable No.-> %d\n",i);
      
      fprintf(m_rep_ptr,"No. of bits assigned to it ->%d\n", m_vlen[i]);
      
      fprintf(m_rep_ptr,"Lower limits on %dth variable-> %f\n",i, m_lim_b[i][0]);
      
      fprintf(m_rep_ptr,"Upper limits on %dth variable ->%f\n",i, m_lim_b[i][1]);
   }
  
   for(i = 0; i < m_nvar; i++)
   {      
      fprintf(m_rep_ptr,"Real-coded variable No.-> %d\n",i);
      
      fprintf(m_rep_ptr,"Lower limits on %dth variable-> %f\n",i,m_lim_r[i][0]);
      
      fprintf(m_rep_ptr,"Upper limits on %dth variable ->%f\n",i,m_lim_r[i][1]);
      if (m_ans == 1)
      {
         fprintf(m_rep_ptr,"Variable bounds are rigid\n");
      }
      else 
      {
         fprintf(m_rep_ptr,"Variable bounds are not rigid\n");
      }
   }/* end for() */

   if (m_nchrom > 0)
   {  
      if(m_optype == 1)
      {
         fprintf(m_rep_ptr,"X-over on binary string is SINGLE POINT X-OVER\n");
      }
      if (m_optype == 2)
      {
         fprintf(m_rep_ptr,"X-over on binary strings is UNIFORM X-OVER \n");
      }
   }/* end if() */
   if (m_nvar > 0) 
   {
      fprintf(m_rep_ptr,"Crossover parameter in the SBX operator is %f\n",m_di);
   }
   fprintf(m_rep_ptr,"Cross-over Probability ->%f\n",m_pcross);
  
   if (m_nchrom > 0)
   {
      fprintf(m_rep_ptr,"Mutation Probability for binary strings -> %f\n", m_pmut_b);
   }

   if (m_nvar > 0)
   {
      fprintf(m_rep_ptr,"Mutation Probability for real-coded vectors -> %f\n", m_pmut_r);
   }
   fprintf(m_rep_ptr,"Random Seed ->%f\n", m_seed);
} /* end InitFromFile() */

/******************************************************************************
crossover()
******************************************************************************/
void NSGAII::crossover(NSGAII_Population *new_pop_ptr,NSGAII_Population *mate_pop_ptr)
{
   int i,k,n,y,mating_site,*par1,*par2,*chld1,*chld2,c;
   float rnd;

   rnd = randomperc();  

   new_pop_ptr->ind_ptr=&(new_pop_ptr->ind[0]);
  
   mate_pop_ptr->ind_ptr=&(mate_pop_ptr->ind[0]); 

   for(i = 0,y = 0,n = 0; i < m_popsize/2; i++)
   {
      new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[n]);
      chld1=&(new_pop_ptr->ind_ptr->genes[0]);
      n = n+1;
      
      new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[n]);
      chld2=&(new_pop_ptr->ind_ptr->genes[0]);
      n = n+1;
      
      mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[y]);
      par1 = &(mate_pop_ptr->ind_ptr->genes[0]);
      y = y+1;
      
      mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[y]);
      par2 = &(mate_pop_ptr->ind_ptr->genes[0]); 
      y = y+1;
      
      rnd = randomperc();
      if (rnd < m_pcross)
      {
         m_ncross++;
         rnd = randomperc();
         c = (int)floor(rnd*(m_chrom+10));
         mating_site = c;

         if(mating_site >= m_chrom)
         {
            mating_site = mating_site/2;
         }
	  
         for(k = 0;k < m_chrom; k++)
         {
            if(k > mating_site-1)
            {
               *chld1++ = *par2++;
               *chld2++ = *par1++;
            }
            else
            {
               *chld1++ = *par1++;
               *chld2++ = *par2++;
            }
         }/* end for(k) */
      }/* end if() */
      else 
      {
         for (k = 0; k < m_chrom; k++)
         {
            *chld1++ = *par1++;
            *chld2++ = *par2++;
         }
      }/* end else() */
   }/* end for(i, y, n) */
}/* end crossover() */

/******************************************************************************
decode()

decode a chromosome to get real values
******************************************************************************/
void NSGAII::decode(NSGAII_Population *pop_ptr)
{
   float *real_ptr;
   int i,sum,b,k,c,d,*gene_ptr,m,x;

   pop_ptr->ind_ptr = &(pop_ptr->ind[0]);
  
   for(i = 0; i < m_popsize; i++)
   {
      real_ptr = &(pop_ptr->ind_ptr->xbin[0]);
      gene_ptr = &(pop_ptr->ind_ptr->genes[0]);

      for(m = 0; m < m_nchrom; m++)
      {
         /*-------------------------------------------------------
         Finding out the co-efficient 2 to the power of 
         (l-1) where l is the no of bits assigned to this variable
	    
         For More Info Study DEB's Book
         -------------------------------------------------------*/
         sum = 0;
         for(k = 0; k < m_vlen[m]; k++)
         {
            b = *gene_ptr;
            d = m_vlen[m] - k - 1;
            c = (int)pow((float)2,d);
            sum = sum + c * b;
            gene_ptr++;
         }/* end for(k) */
	  
         x = m_vlen[m];
         m_coef[m] = pow((float)2,x) - 1;
         *real_ptr =(float)(m_lim_b[m][0] + (sum/m_coef[m])*(m_lim_b[m][1] - m_lim_b[m][0]));
         real_ptr++;
      }/* end for(m) */
      pop_ptr->ind_ptr = &(pop_ptr->ind[i+1]);
   }/* end for(i) */
}/* end decode() */

/******************************************************************************
dumfitness()

Assigns the dummyfitness to the individuals at different ranks
******************************************************************************/
void NSGAII::dumfitness(NSGAII_Population *pop_ptr)
{
   int **ptr; /* [maxpop] */
   int i,j,k,m,rnk1;
  
   FILE *fptr;
   fptr = fopen("/tmp/Dft.dat","w");
   fprintf(fptr," S no F1\tF2\tRank\n");

   ptr = (int **)malloc(sizeof(int *) * m_maxpop);
  
   for(i = 0; i < m_popsize; i++)
   {
      fprintf(fptr,"%d  %f %f %d\n",i+1,pop_ptr->ind[i].fitness[0],pop_ptr->ind[i].fitness[1],pop_ptr->ind[i].rank);
   }

   m = pop_ptr->maxrank;
  
   for(i = 0; i < m; i++)
   {
      m_r_arr[i] = pop_ptr->rankno[i]; /*Put the initial individuals at this rank to zero*/
      fprintf(fptr,"Rank = %d No of Ind = %d\n",i+1,m_r_arr[i]);
   }
  
   for(i = 0; i < m_popsize; i++)
   {   
      pop_ptr->ind[i].cub_len = 0; /*Initialise the dummy fitness to zero*/
   }
      
   for(i = 0; i < m; i++)
   {
      ptr[i] = &(m_rec_arr[i][0]);
   }
  
   for(i = 0; i < m_popsize; i++)
   {
      rnk1 = pop_ptr->ind[i].rank;
      *ptr[rnk1-1]++ = i;
   }
  
   fprintf(fptr,"The rankarray\n");
   for(i = 0 ; i < m; i++)
   {
      fprintf(fptr,"Rank = %d\t",i+1);
      k = m_r_arr[i];
      for(j = 0; j < k; j++)
      {
         fprintf(fptr,"%d ", m_rec_arr[i][j]);
      } 
      fprintf(fptr,"\n");
   }
  
   for(i = 0; i < m; i++)
   {
      share1(pop_ptr,i+1);
   }
  
   fclose(fptr);
   free(ptr);
}/* end dumfitness() */

/******************************************************************************
share1()
******************************************************************************/
void NSGAII::share1(NSGAII_Population *pop_ptr,int rnk)
{
   float ** length; /* [maxpop][2] */
   float max;
   int i,j,m1,a ;

   length = (float **)malloc(sizeof(float *)*m_maxpop);
   for(i = 0; i < m_maxpop; i++)
   {
      length[i] = (float *)malloc(sizeof(float)*2);
   }  

   m1 = m_r_arr[rnk-1];
  
   for(j = 0; j < m_nfunc; j++)
   {
      for(i = 0; i < m1; i++)
      {
         m_fpara[i][0] = 0;
         m_fpara[i][1] = 1;
      }/* end for(i) */
      
      for(i = 0; i < m1; i++)
      {
         a = m_rec_arr[rnk-1][i];
         m_fpara[i][0] = (float)a;
         m_fpara[i][1] = pop_ptr->ind[a].fitness[j];
      }/* end for(i) */
      
      sort1(m1);  /*Sort the array of fitnesses(in ascending order)*/
      
      max = m_fpara[m1-1][1];   /*Find the maximum length*/
      
      for(i = 0; i < m1; i++)
      {
         if(i == 0 || i == (m1-1))
         { 
            /*Put the ends fitness more to preserve them*/
            length[i][0] = m_fpara[i][0];
            length[i][1] = 100*max;
         }/* end if() */
         else
         {
            /*Find the dimension of the cuboid*/
            length[i][0] = m_fpara[i][0];
            length[i][1] = (float)fabs(m_fpara[i+1][1] - m_fpara[i-1][1]);
         }/* end else() */
      }/* end for(i) */
      
      for(i = 0; i < m1; i++)
      {
         a = (int)length[i][0];
         pop_ptr->ind[a].cub_len += length[i][1];
      }/* end for(i) */
   }/* end for(j) */

   for(i = 0; i < m_maxpop; i++)
   {
      free(length[i]);
   }  
   free(length);
}/* end share1() */
 
/******************************************************************************
sort1()

Sorting the length in ascending order
******************************************************************************/
void NSGAII::sort1(int m1)
{
   float temp,temp1; 
   int i1,k1;
   
   for(k1 = 0; k1 < (m1 - 1); k1++)
   {
      for(i1 = (k1 + 1); i1 < m1; i1++)
      {
         if(m_fpara[k1][1] > m_fpara[i1][1])
         {
            temp = m_fpara[k1][1];
            temp1 = m_fpara[k1][0];
            m_fpara[k1][1] = m_fpara[i1][1];
            m_fpara[k1][0] = m_fpara[i1][0];
            m_fpara[i1][1] = temp;
            m_fpara[i1][0] = temp1;
         }/* end if(fpara) */
      }/* end for(i1) */
   }/* end for(k1) */
}/* end sort1() */

/******************************************************************************
func()

Evaluate the objective functions
******************************************************************************/
void NSGAII::func(NSGAII_Population *pop_ptr)
{ 
   float * realx_ptr; /* Pointer to the array of x values */
   float * binx_ptr; /* Pointer to the binary variables */
   float * fitn_ptr; /* Pointer to the array of fitness function */
   double * x; /* [2*maxvar] problem variables */
   double * f; /* [maxfun] array of fitness values */
   float * err_ptr; /* Pointer to the error */
   float * cstr; /* [maxcons] */
   int result;
   int fi;
   double * Xost;
   double * Fost;

   int i,j,k; 
   float error, cc;

   x =(double *)malloc(sizeof(double) * 2 * m_maxvar);
   f = (double *)malloc(sizeof(double) * m_maxfun);
   cstr = (float *)malloc(sizeof(float) * m_maxcons);

   pop_ptr->ind_ptr = &(pop_ptr->ind[0]);

   /*Initializing the max rank to zero*/
   pop_ptr->maxrank = 0;
   for(i = 0 ; i < m_popsize; i++)
   {
      Xost = new double[m_maxvar];
      Fost = new double[m_maxfun];

      pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
      realx_ptr = &(pop_ptr->ind_ptr->xreal[0]);
      binx_ptr = &(pop_ptr->ind_ptr->xbin[0]);
      
      // Real-coded variables
      for(j = 0; j < m_nvar; j++)
      {  
         x[j] = (double)(*realx_ptr++);
         Xost[j] = x[j];
      }

      // Binary-codced variables
      for(j = 0; j < m_nchrom; j++)
      { 
         x[m_nvar+j] = (double)(*binx_ptr++);
         Xost[j] = x[m_nvar+j];
      }
      
      fitn_ptr = &(pop_ptr->ind_ptr->fitness[0]);
      err_ptr = &(pop_ptr->ind_ptr->error);

      /*================================================================
      DO NOT CHANGE ANYTHING ABOVE
      ================================================================*/
      
      /*----------------------------------------------------------------
      CODE YOUR OBJECTIVE FUNCTIONS HERE
   
      All functions must be of minimization type, negate maximization 
      functions.

      Start Coding Your Function From This Point
      -----------------------------------------------------------------*/
      m_pModel->GetParamGroupPtr()->WriteParams(Xost);

      m_pModel->Execute(f, m_maxfun);
      for(fi = 0; fi < m_maxfun; fi++)
      {
         Fost[fi] = f[fi];
      }
      result = UpdateLists(Xost, m_maxvar, Fost, m_maxfun);
      if (result == ARCHIVE_NON_DOM)
      {
         WriteInnerEval(i + 1, m_popsize, '+');
      }/* end if() */
      else
      {
         WriteInnerEval(i + 1, m_popsize, '-');
      }/* end else() */

      /*=========End Your Coding Upto This Point===============*/

      /******************************************************************
      Put The Constraints Here. For example:
         g(x) >= 0 type (normalize g(x) as in the cstr[1] below)                         
      ******************************************************************/      
      // cstr[0] = (float)(x[0]*x[0]+x[1]*x[1]-1.0-0.1*cos(16.0*atan(x[0]/x[1])));
      // cstr[1] = (float)((-NSGAII_SQUARE(x[0]-0.5) - NSGAII_SQUARE(x[1]-0.5) + 0.5)/0.5);
      
      /*================================================================
      Constraints Are Coded Up To Here. 

      DO NOT CHANGE ANYTHING BELOW 
      =================================================================*/
      for(k = 0; k < m_nfunc; k++)
      {
         *fitn_ptr++ = (float)(f[k]);
      }/* end for(k) */
      
      for (k = 0; k < m_ncons; k++)
      {
         pop_ptr->ind_ptr->constr[k] = cstr[k];
      }/* end for(k) */

      error = 0.0;

      for (k = 0; k < m_ncons; k++)
      {
         cc = cstr[k];
         if(cc < 0.0)
         {
            error = error - cc;
         }
      }/* end for(k) */
      *err_ptr = error;
   }/* for (i) */
  
   /*----------------------------------------------------------------
   RANKING 
   ----------------------------------------------------------------*/
   if(m_ncons == 0)
   {
      ranking(pop_ptr);
   }
   else
   {
      rankcon(pop_ptr);
   }

   free(x);
   free(f);
   free(cstr);
}/* end func() */

/******************************************************************************
func_parallel()

Evaluate the objective functions in parallel
******************************************************************************/
void NSGAII::func_parallel(NSGAII_Population *pop_ptr, int myrank, int nprocs)
{ 
   float * realx_ptr; /* Pointer to the array of x values */
   float * binx_ptr; /* Pointer to the binary variables */
   float * fitn_ptr; /* Pointer to the array of fitness function */
   double * x; /* [2*maxvar] problem variables */
   double * f; /* [maxfun] array of fitness values */
   float * err_ptr; /* Pointer to the error */
   float * cstr; /* [maxcons] */
   int result;
   int fi;
   double * Xost;
   double * Fost;

   int i,j,k; 
   float error, cc;

   x =(double *)malloc(sizeof(double) * 2 * m_maxvar);
   f = (double *)malloc(sizeof(double) * m_maxfun);
   cstr = (float *)malloc(sizeof(float) * m_maxcons);

   pop_ptr->ind_ptr = &(pop_ptr->ind[0]);

   /*Initializing the max rank to zero*/
   pop_ptr->maxrank = 0;
   for(i = 0 ; i < m_popsize; i++)
   {
      Xost = new double[m_maxvar];
      Fost = new double[m_maxfun];

      pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
      realx_ptr = &(pop_ptr->ind_ptr->xreal[0]);
      binx_ptr = &(pop_ptr->ind_ptr->xbin[0]);
      
      // Real-coded variables
      for(j = 0; j < m_nvar; j++)
      {  
         x[j] = (double)(*realx_ptr++);
         Xost[j] = x[j];
      }

      // Binary-codced variables
      for(j = 0; j < m_nchrom; j++)
      { 
         x[m_nvar+j] = (double)(*binx_ptr++);
         Xost[j] = x[m_nvar+j];
      }
      
      fitn_ptr = &(pop_ptr->ind_ptr->fitness[0]);
      err_ptr = &(pop_ptr->ind_ptr->error);

      /*================================================================
      DO NOT CHANGE ANYTHING ABOVE
      ================================================================*/
      
      /*----------------------------------------------------------------
      CODE YOUR OBJECTIVE FUNCTIONS HERE
   
      All functions must be of minimization type, negate maximization 
      functions.

      Start Coding Your Function From This Point
      -----------------------------------------------------------------*/
      m_pModel->GetParamGroupPtr()->WriteParams(Xost);

      m_pModel->Execute(f, m_maxfun);
      for(fi = 0; fi < m_maxfun; fi++)
      {
         Fost[fi] = f[fi];
      }
      result = UpdateLists(Xost, m_maxvar, Fost, m_maxfun);
      if (result == ARCHIVE_NON_DOM)
      {
         WriteInnerEval(i + 1, m_popsize, '+');
      }/* end if() */
      else
      {
         WriteInnerEval(i + 1, m_popsize, '-');
      }/* end else() */

      /*=========End Your Coding Upto This Point===============*/

      /******************************************************************
      Put The Constraints Here. For example:
         g(x) >= 0 type (normalize g(x) as in the cstr[1] below)                         
      ******************************************************************/      
      // cstr[0] = (float)(x[0]*x[0]+x[1]*x[1]-1.0-0.1*cos(16.0*atan(x[0]/x[1])));
      // cstr[1] = (float)((-NSGAII_SQUARE(x[0]-0.5) - NSGAII_SQUARE(x[1]-0.5) + 0.5)/0.5);
      
      /*================================================================
      Constraints Are Coded Up To Here. 

      DO NOT CHANGE ANYTHING BELOW 
      =================================================================*/
      for(k = 0; k < m_nfunc; k++)
      {
         *fitn_ptr++ = (float)(f[k]);
      }/* end for(k) */
      
      for (k = 0; k < m_ncons; k++)
      {
         pop_ptr->ind_ptr->constr[k] = cstr[k];
      }/* end for(k) */

      error = 0.0;

      for (k = 0; k < m_ncons; k++)
      {
         cc = cstr[k];
         if(cc < 0.0)
         {
            error = error - cc;
         }
      }/* end for(k) */
      *err_ptr = error;
   }/* for (i) */
  
   /*----------------------------------------------------------------
   RANKING 
   ----------------------------------------------------------------*/
   if(m_ncons == 0)
   {
      ranking(pop_ptr);
   }
   else
   {
      rankcon(pop_ptr);
   }

   free(x);
   free(f);
   free(cstr);
}/* end func_parallel() */

/******************************************************************************
init()

initialize the population
******************************************************************************/
void NSGAII::init(NSGAII_Population *pop_ptr)
{
   int i,j;
   float d;
   pop_ptr->ind_ptr = &(pop_ptr->ind[0]); 
  
   /*Loop Over the population size*/
   for (i = 0; i < m_popsize; i++)
   { 
      /*Loop over the chromosome length*/
      for (j = 0; j < m_chrom; j++)
      {
         /*-----------------------------------------------
         Generate a Random No. if it is less than 0.5 it 
         generates a 0 in the string otherwise 1
         -----------------------------------------------*/
         d = randomperc();
         if(d >= 0.5)
         {
            pop_ptr->ind_ptr->genes[j] = 1;
         }
         else
         {
            pop_ptr->ind_ptr->genes[j] = 0;
         } 
      }/* end for() */
      pop_ptr->ind_ptr = &(pop_ptr->ind[i+1]);
   }/* end for() */
   pop_ptr->ind_ptr = &(pop_ptr->ind[0]); 
}/* end init() */

/******************************************************************************
keepalive()
******************************************************************************/
void NSGAII::keepalive(NSGAII_Population *pop1_ptr,NSGAII_Population *pop2_ptr,NSGAII_Population *pop3_ptr,int gen)
{
   int i,j,jj,k,m,l,rec;  
   int st,pool,poolf,sel;
   int *gene1_ptr, *gene2_ptr;
  
   float *gene3_ptr,*gene4_ptr,*xbin1_ptr,*xbin2_ptr;

   /* Forming the global mating pool */
   for(i = 0; i < m_popsize; i++)
   {
      if(m_nchrom > 0)
      {
         /*Binary Coded GA genes are copied*/
         for(k = 0; k < m_chrom; k++)
         {
            m_globalpop.genes[i][k]=pop1_ptr->ind[i].genes[k];
            m_globalpop.genes[i+m_popsize][k] = pop2_ptr->ind[i].genes[k];
         }
         for(k = 0; k < m_nchrom; k++)
         {
            m_globalpop.xbin[i][k] = pop1_ptr->ind[i].xbin[k];
            m_globalpop.xbin[i+m_popsize][k] = pop2_ptr->ind[i].xbin[k];
         }
      }/* end if() */
      if(m_nvar > 0)
      {
         /*For Real Coded GA x values are copied */
         for(k = 0; k < m_nvar; k++)
         {
            m_globalpop.xreal[i][k] = pop1_ptr->ind[i].xreal[k];
            m_globalpop.xreal[i+m_popsize][k] = pop2_ptr->ind[i].xreal[k];
         }
      }
      
      /*Fitness is copied to the global pool */
      for(l = 0; l < m_nfunc; l++)
      {
         m_globalpop.fitness[i][l] = pop1_ptr->ind[i].fitness[l];
         m_globalpop.fitness[i+m_popsize][l] = pop2_ptr->ind[i].fitness[l];
      }
       
      /*Initial;ising the dummyfitness to zero */
      m_globalpop.cub_len[i] = 0;
      m_globalpop.cub_len[i+m_popsize] = 0;
      m_globalpop.error[i] = pop1_ptr->ind[i].error;
      m_globalpop.error[i+m_popsize] = pop2_ptr->ind[i].error;

      for (jj=0; jj < m_ncons; jj++)
      {
         m_globalpop.constr[i][jj] = pop1_ptr->ind[i].constr[jj];
         m_globalpop.constr[i+m_popsize][jj] = pop2_ptr->ind[i].constr[jj];
      }
   }/* end for(i) */
  
   m_global_pop_ptr = &(m_globalpop);
  
   /*Finding the global ranks */
   if(m_ncons == 0)
   {
      grank(gen);
   }
   else
   {
      grankc(gen);
   }
   m = m_globalpop.maxrank;
  
   /* Sharing the fitness to get the dummy fitness */
   for(i = 0; i < m; i++)
   {
      gshare(i+1);
   }
  
   poolf = m_popsize;
   pool = 0;
  
  
   /*Initializing the flags of population to zero */
   for(i = 0;i < 2 * m_popsize;i++)
   {
      m_globalpop.flag[i] = 0;
   }
   // decide which all solutions belong to the pop3 
   rec = 0;
   st = 0;
   for(i = 0; i < m; i++)
   {
      /*    Elitism Applied Here     */
      st = pool;
      pool += m_globalpop.rankno[i];
      
      if(pool <= m_popsize)
      {
         for(k = 0;k < 2*m_popsize ;k++)
         {
            if(m_globalpop.rank[k] == i+1)
            {
               m_globalpop.flag[k] = 1;
            }
         }/* end for() */

         pop3_ptr->rankno[i] = m_globalpop.rankno[i];
      }/* end if() */
      else
      {
         sel = m_popsize - st;
         m_Lastrank = i+1;
         pop3_ptr->rankno[i] = sel;
         gsort(i+1,sel);
         break;
      }/* end else() */
   }/* end for(i) */

   k = 0;
   for(i = 0,k = 0; i < (2 * m_popsize) && k < m_popsize; i++)
   {
      if(m_nchrom > 0)
      {
         if(m_globalpop.flag[i] == 1)
         {
            gene1_ptr = &(m_globalpop.genes[i][0]);
            xbin1_ptr = &(m_globalpop.xbin[i][0]);
            pop3_ptr->ind_ptr = &(pop3_ptr->ind[k]);
            gene2_ptr = &(pop3_ptr->ind_ptr->genes[0]);
            xbin2_ptr = &(pop3_ptr->ind_ptr->xbin[0]);
	      
            for(j = 0 ; j < m_chrom; j++)
            {
               *gene2_ptr++ = *gene1_ptr++;
            }
            for (j=0; j < m_nchrom; j++)
            { 
               *xbin2_ptr++ = *xbin1_ptr++;
            }
         }/* end if() */
      }/* end if() */
      if (m_nvar > 0)
      {
         if(m_globalpop.flag[i] == 1)
         {
            gene3_ptr = &(m_globalpop.xreal[i][0]);
            pop3_ptr->ind_ptr = &(pop3_ptr->ind[k]);
            gene4_ptr = &(pop3_ptr->ind_ptr->xreal[0]);

            for(j = 0; j < m_nvar; j++)
            {
               *gene4_ptr++ = *gene3_ptr++;
            } 
         }/* end if() */
      }/* end if() */
      if(m_globalpop.flag[i] == 1)
      {
         for(j = 0; j < m_nfunc; j++)
         {
            pop3_ptr->ind[k].fitness[j] = m_globalpop.fitness[i][j];
         }
         pop3_ptr->ind[k].cub_len = m_globalpop.cub_len[i];
         if(m_ncons != 0)
         {
            pop3_ptr->ind[k].error = m_globalpop.error[i];
         }
         for (jj=0; jj < m_ncons; jj++)
         { 
            pop3_ptr->ind[k].constr[jj] = m_globalpop.constr[i][jj];
         }
         pop3_ptr->ind[k].rank = m_globalpop.rank[i];
         k++;  // increment the pop3 counter
      }/* end if() */
   }/* end for(i) */
  
   pop3_ptr->maxrank = m_Lastrank;
}/* end keepalive() */

/******************************************************************************
grank()
******************************************************************************/
void NSGAII::grank(int gen)
{
   int i,j,k,rnk,val,nondom,popsize1,q;
   int * gflg; /* [2*maxpop] */
   float *ptr1,*ptr2;
   FILE *gr;
   gr = fopen("NSGAII_g_rank_record.out","a");
   fprintf(gr,"Genration no. = %d\n",gen);
   /*----------------------------* RANKING *---------------------------------*/
   rnk = 0;
   nondom = 0;
   popsize1 = 2*m_popsize;

   gflg = (int *)malloc(sizeof(int)*2*m_maxpop);
  
   for(i = 0;i < popsize1;i++)
   {
      gflg[i] = 2;
   }
  
   for(k = 0;k < popsize1;k++)
   {
      q =  0;
      for(j = 0; j < popsize1; j++)
      {
         if (gflg[j] != 1) 
         { 
            break;
         }
      }
      if(j == popsize1) 
      {
         break;
      }
      rnk = rnk +1;
      for(j = 0; j < popsize1; j++)
      {
         if(gflg[j] == 0)
         { 
            gflg[j] = 2;
         }
      }
      for(i = 0; i < popsize1; i++)
      {
         if(gflg[i] != 1 && gflg[i] != 0) 
         {
            ptr1 = &(m_global_pop_ptr->fitness[i][0]);
            for(j = 0; j < popsize1; j++)
            {
               if( i!= j)
               {
                  if(gflg[j] != 1)
                  {
                     ptr2 = &(m_global_pop_ptr->fitness[j][0]);
                     val = indcmp1(ptr1,ptr2);
                     if( val == 2)
                     { 
                        gflg[i] = 0;/* individual 1 is dominated */
                        break;
                     }
                     if(val == 1)
                     {
                        gflg[j] = 0;/* individual 2 is dominated */
                     }
                     if(val == 3)
                     {
                        nondom++;/* individual 1 & 2 are non dominated */
                        if(gflg[j] != 0)
                        {
                           gflg[j] = 3;
                        }
                     }
                  }
               }
            }
            if( j == popsize1)
            {
               m_global_pop_ptr->rank[i] = rnk;
               gflg[i] = 1;
               m_global_pop_ptr->rankar[rnk-1][q] =  i;
               q++;
            }
         }
      }
      m_global_pop_ptr->rankno[rnk-1] = q;
   } 
   m_global_pop_ptr->maxrank = rnk;

   fprintf(gr,"   RANK     No Of Individuals\n");
   for(i = 0; i < rnk; i++)
   {
      fprintf(gr,"\t%d\t%d\n",i+1,m_globalpop.rankno[i]);
   }

   fclose(gr);
   free(gflg);
}/* end grank() */

/******************************************************************************
grankc()
******************************************************************************/
void NSGAII::grankc(int gen)
{
   int i,j,k,rnk,val,nondom,popsize1,q;
   int * gflg; /* [2*maxpop] */
   float *ptr1,*ptr2;
   float *err_ptr1,*err_ptr2;
   FILE *gr;
   gr = fopen("NSGAII_g_rank_record.out","a");
   fprintf(gr,"Genration no. = %d\n",gen);

   gflg = (int *)malloc(sizeof(int)*2*m_maxpop);

   /*----------------------------* RANKING *---------------------------------*/
   rnk = 0;
   nondom = 0;
   popsize1 = 2*m_popsize;
   m_min_fit = (float)popsize1;
   m_delta_fit = (float)(0.1 *popsize1);
   for(i = 0; i < popsize1; i++)
   { 
      gflg[i] = 2;
   }
   for(k = 0; k < popsize1; k++)
   {
      q =  0;
      for(j = 0; j < popsize1; j++)
      {
         if (gflg[j] != 1) 
         {
            break;
         }
      }
      if(j == popsize1) 
      {
         break;
      }
      rnk = rnk +1;
      for(j = 0; j < popsize1; j++)
      {
         if(gflg[j] == 0) 
         { 
            gflg[j] = 2;
         }
      }
      for(i = 0; i < popsize1; i++)
      {
         if(gflg[i] != 1 && gflg[i] != 0) 
         {
            ptr1 = &(m_global_pop_ptr->fitness[i][0]);
            err_ptr1 = &(m_global_pop_ptr->error[i]);
            for(j = 0; j < popsize1; j++)
            {
               if( i!= j)
               {
                  if(gflg[j] != 1)
                  {
                     ptr2 = &(m_global_pop_ptr->fitness[j][0]);
                     err_ptr2 = &(m_global_pop_ptr->error[j]);

                     if(*err_ptr1 < 1.0e-6 && *err_ptr2 > 1.0e-6)
                     {
                        /* first feasible second individaul is infeasible*/
                        gflg[j] = 0;
                     }
                     else
                     {
                        if(*err_ptr1 >1.0e-6 && *err_ptr2 < 1.0e-6)
                        {
                           /*first individual is infeasible and second is feasible*/
                           gflg[i] = 0;
                           break;
                        }
                        else /*both feasible or both infeasible*/
                        {
                           if(*err_ptr1 > *err_ptr2)
                           {
                              gflg[i] = 0;
                              /*first individual is more infeasible*/
                              break;
                           }
                           else
                           {
                              if(*err_ptr1 < *err_ptr2)
                              { 
                                 gflg[j] = 0;
                              }
                              /*second individual is more infeasible*/
                              else
                              {
                                 val = indcmp1(ptr1,ptr2);
                                 if( val == 2)
                                 { 
                                    gflg[i] = 0;
                                    /* individual 1 is dominated */
                                    break;
                                 }
                                 if(val == 1)
                                 {
                                    gflg[j] = 0;
                                    /* individual 2 is dominated */
                                 }
                                 if(val == 3)
                                 {
                                    nondom++;/* individual 1 & 2 are non dominated */
                                    if(gflg[j] != 0) 
                                    { 
                                       gflg[j] = 3;
                                    }
                                 }/* end if(val) */
                              }/* end else() */
                           }/* end else() */
                        }/* ens else() */
                     }/* end else() */
                  }/* end if(gflg[j] != 1) */
               }/* end if(i != j)*/
            }/* end for(j) */
            if(j == popsize1)
            {
               m_global_pop_ptr->rank[i] = rnk;
               gflg[i] = 1;
               m_global_pop_ptr->rankar[rnk-1][q] =  i;
               q++;
            }/* end if() */
         }/* end if(flags) */
      }/* end for(i) */
      m_global_pop_ptr->rankno[rnk-1] = q;
   }/* end for(k) */ 
   m_global_pop_ptr->maxrank = rnk;
   fprintf(gr,"   RANK     No Of Individuals\n");

   for(i = 0;i < rnk;i++)
   {
      fprintf(gr,"\t%d\t%d\n",i+1,m_globalpop.rankno[i]);
   }

   fclose(gr);
   free(gflg);
}/* end grankc() */

/******************************************************************************
Create_NSGAII_GlobPop()

Create global pop data struct.
******************************************************************************/
void NSGAII::Create_NSGAII_GlobPop(NSGAII_GlobPop * pGlobalPop)
{
   int i;

   pGlobalPop->rankar = new int *[2*m_maxpop];
   for(i = 0; i < 2*m_maxpop; i++)
   {
      pGlobalPop->rankar[i] = new int[2*m_maxpop];
   }
   pGlobalPop->rankno = new int[2*m_maxpop];
   pGlobalPop->rank = new int[2*m_maxpop];
   pGlobalPop->flag = new int[2*m_maxpop];

   pGlobalPop->genes = new int *[2*m_maxpop];
   for(i = 0; i < 2*m_maxpop; i++)
   {
      pGlobalPop->genes[i] = new int[m_maxchrom];
   }

   pGlobalPop->fitness = new float *[2*m_maxpop];
   for(i = 0; i < 2*m_maxpop; i++)
   {
      pGlobalPop->fitness[i] = new float[m_maxfun];
   }

   pGlobalPop->cub_len = new float[2*m_maxpop];

   pGlobalPop->xreal = new float *[2*m_maxpop];
   for(i = 0; i < 2*m_maxpop; i++)
   {
      pGlobalPop->xreal[i] = new float[m_maxvar];
   }

   pGlobalPop->xbin = new float *[2*m_maxpop];
   for(i = 0; i < 2*m_maxpop; i++)
   {
      pGlobalPop->xbin[i] = new float[m_maxvar];
   }

   pGlobalPop->error = new float[2*m_maxpop];

   pGlobalPop->constr = new float *[2*m_maxpop];
   for(i = 0; i < 2*m_maxpop; i++)
   {
      pGlobalPop->constr[i] = new float[m_maxcons];
   }
}/* end Create_NSGAII_GlobPop() */

/******************************************************************************
Create_NSGAII_Population()

Create population data struct.
******************************************************************************/
void NSGAII::Create_NSGAII_Population(NSGAII_Population * pPop)
{
   int i;
   pPop->rankrat = new float[m_maxpop]; /* [maxpop] Rank Ratio */
   pPop->rankno = new int[m_maxpop]; /* [maxpop] Individual at different ranks */
   pPop->ind = new NSGAII_Individual[m_maxpop]; /* [maxpop] Different Individuals */
   for(i = 0; i < m_maxpop; i++)
   {
      Create_NSGAII_Individual(&(pPop->ind[i]));
   }/* end for() */
}/* end Create_NSGAII_Population() */

/******************************************************************************
Create_NSGAII_Individual()

Create individuak data struct.
******************************************************************************/
void NSGAII::Create_NSGAII_Individual(NSGAII_Individual * pInd)
{
   pInd->genes = new int[m_maxchrom]; /* [maxchrom] bianry chromosome */
   pInd->gener = new float[m_maxchrom]; /* [maxchrom] real chromosomes */
   pInd->xreal = new float[m_maxvar]; /* [maxvar] list of real variables */
   pInd->xbin = new float[m_maxvar]; /* [maxvar] list of decoded value of the chromosome */
   pInd->fitness = new float[m_maxfun]; /* [maxfun] Fitness values */
   pInd->constr = new float[m_maxcons]; /* [maxcons] Constraints values */
}/* end Create_NSGAII_Individual() */

/******************************************************************************
Delete_NSGAII_GlobPop()

Cleanup global pop data struct.
******************************************************************************/
void NSGAII::Delete_NSGAII_GlobPop(NSGAII_GlobPop * pGlobalPop)
{
   int i;

   for(i = 0; i < 2*m_maxpop; i++)
   {
      delete [] pGlobalPop->rankar[i];
      delete [] pGlobalPop->genes[i];
      delete [] pGlobalPop->fitness[i]; /* [2*maxpop][maxfun] Fitness function values for the different individuals */
      delete [] pGlobalPop->xreal[i]; /* [2*maxpop][maxvar] value of the decoded variables for different individuals */
      delete [] pGlobalPop->xbin[i]; /* [2*maxpop][maxvar] binray-coded variables */
      delete [] pGlobalPop->constr[i]; /* [2*maxpop][maxcons] */
   }
   delete [] pGlobalPop->rankar;
   delete [] pGlobalPop->genes;
   delete [] pGlobalPop->fitness; /* [2*maxpop][maxfun] Fitness function values for the different individuals */
   delete [] pGlobalPop->xreal; /* [2*maxpop][maxvar] value of the decoded variables for different individuals */
   delete [] pGlobalPop->xbin; /* [2*maxpop][maxvar] binray-coded variables */
   delete [] pGlobalPop->constr; /* [2*maxpop][maxcons] */

   delete [] pGlobalPop->rankno; /* [2*maxpop] record of no. of individuals at a particular rank */     
   delete [] pGlobalPop->rank; /* [2*maxpop] rank of different individuals */
   delete [] pGlobalPop->flag; /* [2*maxpop] Setting the flag */

   delete [] pGlobalPop->cub_len; /* [2*maxpop] Dummyfitness */
   delete [] pGlobalPop->error; /* [2*maxpop] Error Values of the individuals */
}/* end Delete_NSGAII_GlobPop() */

/******************************************************************************
indcmp1()
******************************************************************************/
int NSGAII::indcmp1(float *ptr1,float *ptr2)
{
   float * fit1; /* [maxfun] */
   float * fit2; /* [maxfun] */
   int i,value,m,n;

   fit1 = (float *)malloc(sizeof(float)*m_maxfun);
   fit2 = (float *)malloc(sizeof(float)*m_maxfun);

   for(i = 0; i < m_nfunc; i++)
   {
      fit1[i] = *ptr1++;
      fit2[i] = *ptr2++;
   }
   m = 0;
   n = 0;
   while(m < m_nfunc && fit1[m] <= fit2[m]) 
   {
      if((fit2[m] - fit1[m]) < 1e-7)
      {
         n++;
      }
      m++;
   }
   if(m == m_nfunc) 
   {
      if(n == m_nfunc) 
      {
         value = 3;
      }
      else 
      { 
         value = 1; /*value = 1 for dominating*/
      }
   }
   else 
   {
      m = 0;
      n = 0;
      while(m < m_nfunc && fit1[m] >= fit2[m]) 
      {
         if((fit1[m] - fit2[m]) < 1e-7) 
         {
            n++;
         }
         m++;
      }
      if(m == m_nfunc)
      {
         if(n != m_nfunc)
         {
            value = 2; /*value =  2 for dominated */
         }
         else 
         { 
            value = 3;
         }
      }
      else  
      { 
         value = 3; /*value = 3 for incomparable*/
      }   
   }
   free(fit1);
   free(fit2);
   return value;
}/* end indcmp1() */

/******************************************************************************
gsort()

sort the dummyfitness arrays
******************************************************************************/
void NSGAII::gsort(int rnk,int sel)
{
   int i,j,a,q;
   float ** array; /* [2*maxpop][2] */
   float temp,temp1;
  
   array = (float **)malloc(sizeof(float *)*2*m_maxpop);
   for(i = 0; i < 2*m_maxpop; i++)
   {
      array[i] = (float *)malloc(sizeof(float)*2);
   }

   q = m_globalpop.rankno[rnk-1];
  
   for(i = 0; i < q; i++)
   {
      array[i][0] = (float)m_globalpop.rankar[rnk-1][i];
      a = m_globalpop.rankar[rnk-1][i];
      array[i][1] = m_globalpop.cub_len[a];
   }
   for(i = 0; i < q; i++)
   {
      for(j = i+1; j < q; j++)
      {
         if(array[i][1] < array[j][1])
         {
            temp = array[i][1];
            temp1 = array[i][0];
            array[i][1] = array[j][1];
            array[i][0] = array[j][0];
	      
            array[j][1] = temp;
            array[j][0] = temp1;
         }
      }
   }
  
   for(i = 0; i < sel; i++)
   {
      a = (int)array[i][0];
      m_globalpop.flag[a] = 1;
   }

   for(i = 0; i < 2*m_maxpop; i++)
   {
      free(array[i]);
   }
   free(array);
}/* end gsort() */

/******************************************************************************
gshare()
******************************************************************************/
void NSGAII::gshare(int rnk)
{
   float ** length; /* [2*maxpop][2] */
   float max;
   int i,j,m1,a ;

   length = (float **)malloc(sizeof(float *) * (m_maxpop * 2));
   for(i = 0; i < (2*m_maxpop); i++)
   {
      length[i] = (float *)malloc(sizeof(float)*2);
   }  

   m1 = m_globalpop.rankno[rnk-1];
  
   for(j = 0; j < m_nfunc; j++)
   {
      for(i = 0; i < m1; i++)
      {
         m_fpara1[i][0] = 0;
         m_fpara1[i][1] = 0;
      }
      
      for(i = 0; i < m1; i++)
      {
         a = m_globalpop.rankar[rnk-1][i];
         m_fpara1[i][0] = (float)a ;
         m_fpara1[i][1] = m_globalpop.fitness[a][j];
      }
      
      sort(m1); /*Sort the arrays in ascending order of the fitness*/
      
      max = m_fpara1[m1-1][1];
      for(i = 0; i < m1; i++)
      {
         if(i == 0 ||i == (m1-1))
         { 
            length[i][0] = m_fpara1[i][0];
            length[i][1] = 100 * max;
         }
         else
         {
            length[i][0] = m_fpara1[i][0];
            length[i][1] = (float)fabs(m_fpara1[i+1][1] - m_fpara1[i-1][1]);
         }
      }
      for(i = 0; i < m1; i++)
      {
         a = (int)length[i][0];
         m_globalpop.cub_len[a] += length[i][1];
      }
   }

   for(i = 0; i < 2*m_maxpop; i++)
   {
      free(length[i]);
   }  
   free(length);
}/* end gshare() */

/******************************************************************************
sort()
******************************************************************************/
void NSGAII::sort(int m1)
{
   float temp,temp1; 
   int i1,k1;

   for(k1 = 0;k1 < m1-1;k1++)
   {
      for(i1 = k1+1;i1 < m1;i1++)
      {
         if(m_fpara1[k1][1] > m_fpara1[i1][1])
         {
            temp = m_fpara1[k1][1];
            temp1 = m_fpara1[k1][0];
            m_fpara1[k1][1] = m_fpara1[i1][1];
            m_fpara1[k1][0] = m_fpara1[i1][0];
            m_fpara1[i1][1] = temp;
            m_fpara1[i1][0] = temp1;
         }
      }
   }
}/* end sort() */

/******************************************************************************
mutate()

ormulate the mutation routine
******************************************************************************/
void NSGAII::mutate(NSGAII_Population *new_pop_ptr)
{
   int i,*ptr,j;
   float rand1;

   rand1 = randomperc();
   new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);
  
   for(j = 0; j < m_popsize; j++)
   {
      ptr = &(new_pop_ptr->ind_ptr->genes[0]);
      new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j+1]);     
      
      /*Select bit */
      for (i = 0; i < m_chrom; i++)
      {
         rand1 = randomperc();
	  
         /*Check whether to do mutation or not*/
         if(rand1 <= m_pmut_b)
         {
            if(*ptr == 0)
            {
               *ptr = 1;
            }
            else
            {
               *ptr = 0;
            }
            m_nmut++;
         }/* end if() */
         ptr++;
      }/* end for() */
   }/* end for() */
}/* end mutate() */

/******************************************************************************
rankcon()

This also demarkates the different Pareto Fronts
*****************************************************************************/
void NSGAII::rankcon(NSGAII_Population *pop_ptr)
{
   int i,j,k; /* counters */
   int rnk; /* rank */
   int val; /* value obtained after comparing two individuals */
   int nondom; /* no of non dominated members */
   int maxrank1; /* Max rank of the population */
   int * rankarr; /* [maxpop] Array storing the individual number at a rank */
   int q;

   float *ptr1,*ptr2,*err_ptr1,*err_ptr2;

   rankarr = (int *)malloc(sizeof(int) *  m_maxpop);

   /*----------------------------------------------------------------
   RANKING 
   ----------------------------------------------------------------*/

   /*Initializing the ranks to zero*/
   rnk = 0 ;

   nondom = 0 ;
   maxrank1 = 0;

   /*----------------------------------------------------------------
   min_fit is initialize to start distributing the dummy fitness = 
   popsize to the rank one individuals and keeping the record such 
   that the minimum fitness of the better rank individual is always 
   greater than max fitness of the relatively worse rank
   ----------------------------------------------------------------*/
   m_min_fit = (float)m_popsize;


   /*----------------------------------------------------------------
   Difference in the fitness of minimum dummy fitness of better rank 
   and max fitness of the next ranked individuals
   ----------------------------------------------------------------*/
   m_delta_fit = (float)(0.1 * m_popsize);

   /* Initializing all the flags to 2 */
   for(j = 0; j < m_popsize; j++)
   {
      pop_ptr->ind[j].flag = 2;
   }

   q = 0;

   for(k = 0; k < m_popsize; k++,q=0)
   {
      for(j = 0; j < m_popsize; j++)
      {
         if (pop_ptr->ind[j].flag != 1)
         { 
            break; /* Break if all the individuals are assigned a rank */
         }
      }/* end for() */
      if(j == m_popsize)
      {
         break;
      }/* end if() */

      rnk = rnk + 1;

      for( j = 0 ;j < m_popsize; j++)
      {
         if(pop_ptr->ind[j].flag == 0) 
         {
            pop_ptr->ind[j].flag = 2; /* Set the flag of dominated individuals to 2 */
         }
      }/* end for() */

      for(i = 0; i < m_popsize; i++)
      {
         /*Select an individual which rank to be assigned*/

         pop_ptr->ind_ptr = &(pop_ptr->ind[i]);

         if(pop_ptr->ind_ptr->flag != 1 && pop_ptr->ind_ptr->flag != 0) 
         {
            ptr1 = &(pop_ptr->ind_ptr->fitness[0]);
            err_ptr1 = &(pop_ptr->ind_ptr->error);

            for(j = 0;j < m_popsize ; j++)
            {
               /*Select the other individual which has not got a rank*/
               if(i != j)
               {
                  if(pop_ptr->ind[j].flag != 1)
                  {
                     pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
                     ptr2 = &(pop_ptr->ind_ptr->fitness[0]);
                     err_ptr2 = &(pop_ptr->ind_ptr->error);
			  
                     if(*err_ptr1 < 1.0e-6 && *err_ptr2 > 1.0e-6)
                     {
                        /*first ind is feasible second individual is infeasible*/
                        pop_ptr->ind[j].flag = 0;
                     }/* end if() */
                     else
                     {
                        if(*err_ptr1 > 1.0e-6 && *err_ptr2 < 1.0e-6)
                        {
                           /*first individual is infeasible and second is feasible*/
                           pop_ptr->ind[i].flag = 0;
                           break;
                        }/* end if() */
                        else
                        {
                           /*both are feasible or both are infeasible*/
                           if(*err_ptr1 > *err_ptr2)
                           {
                              pop_ptr->ind[i].flag = 0;
                              /*first individual is more infeasible*/
                              break;
                           }/* end if() */
                           else
                           {
                              if(*err_ptr1 < *err_ptr2)
                              {
                                 pop_ptr->ind[j].flag = 0;
                                 /*second individual is more infeasible*/
                              }
                              else
                              {
                                 /*Compare the two individuals for fitness*/
                                 val = indcmp3(ptr1,ptr2);
					  
                                 /*VAL = 2 for dominated individual which rank to be given*/
                                 /*VAL = 1 for dominating individual which rank to be given*/
                                 /*VAL = 3 for non comparable individuals*/
                                 if( val == 2)
                                 { 
                                    pop_ptr->ind[i].flag = 0;
                                    /* individual 1 is dominated */
                                    break;
                                 }
                                 if(val == 1)
                                 {
                                    pop_ptr->ind[j].flag = 0;
                                    /* individual 2 is dominated */
                                 }
                                 if(val == 3)
                                 {
                                    nondom++;
                                    /* individual 1 & 2 are non dominated */
                                    if(pop_ptr->ind[j].flag != 0)
                                    {
                                       pop_ptr->ind[j].flag = 3;
                                    }
                                 }/* end if() */
                              } /* end else() */
                           } /* end else */
                        }/* end else() */
                     }/* end else() */
                  }/* end if() */
               }/* end if(i != j) */
            } /* end for(j) */
            if(j == m_popsize)
            {
               /*Assign the rank and set the flag*/
               pop_ptr->ind[i].rank = rnk;
               pop_ptr->ind[i].flag = 1;
               rankarr[q] = i;
               q++;
            }/* end if() */
         }/* end if(flag) */
      } /* end for(i) */
      pop_ptr->rankno[rnk-1] = q ;
   }/* end for(k) */
   maxrank1 = rnk;
  
   /* Find Max Rank of the population */
   for(i = 0; i < m_popsize; i++)
   {
      rnk = pop_ptr->ind[i].rank;

      if(rnk > maxrank1) 
      { 
         maxrank1 = rnk;
      }
   }/* end for() */

   pop_ptr->maxrank = maxrank1;

   free(rankarr);
}/* end rankcon() */

/******************************************************************************
indcmp3()

Routine Comparing the two individuals
*****************************************************************************/
int NSGAII::indcmp3(float *ptr1,float *ptr2)
{
   float * fit1; /* [maxfun] */
   float * fit2; /* [maxfun] */
   int i,value,m,n;

   fit1 = (float *)malloc(sizeof(float)*m_maxfun);
   fit2 = (float *)malloc(sizeof(float)*m_maxfun);

   for(i = 0;i < m_nfunc ;i++)
   {
      fit1[i] = *ptr1++;
      fit2[i] = *ptr2++;
   }
   m = 0;
   n = 0;
   while(m < m_nfunc && fit1[m] <= fit2[m]) 
   {
      if(fit1[m]== fit2[m]) 
      {
         n++;
      }
      m++;
   }/* end while() */
   if(m == m_nfunc) 
   {
      if(n == m_nfunc) 
      { 
         value = 3;
      }
      else 
      { 
         value = 1; /* value = 1 for dominationg */
      }
   }/* end if() */
   else 
   {
      m = 0;
      n = 0;
      while(m < m_nfunc && fit1[m] >= fit2[m])
      {
         if(fit1[m]== fit2[m]) 
         {
            n++;
         }
         m++;
      }/* end while() */
      if(m == m_nfunc)
      {
         if(n != m_nfunc)
         { 
            value = 2;
         } /* value =  2 for dominated */
         else 
         {
            value = 3;
         }
      }/* end if() */
      else 
      { 
         value = 3; /* value = 3 for incomparable */
      }
   }/* end else() */
  
   free(fit1);
   free(fit2);

   return value;
}/* end indcmp3() */

/******************************************************************************
flip()

Flip a biased coin - true if heads
*****************************************************************************/
int NSGAII::flip(float prob) 
{  
   if(randomperc() <= prob) 
   {
      return(1); 
   }
   else 
   {
      return(0);
   } 
} /* end flip() */
     
/******************************************************************************
noise()

normal noise with specified mean & std dev: mu & sigma
*****************************************************************************/
double NSGAII::noise(double mu, double sigma) 
{
   return((randomnormaldeviate()*sigma) + mu); 
} /* end noise() */
  
/******************************************************************************
randomnormaldeviate()

random normal deviate
*****************************************************************************/
double NSGAII::randomnormaldeviate(void) 
{ 
   return(GaussRandom());
} /* end randomnormaldeviate() */

/******************************************************************************
randomperc()

Fetch a single random number between 0.0 and 1.0 
*****************************************************************************/
float NSGAII::randomperc(void) 
{ 
   return ((float)UniformRandom());
} /* end randomperc() */

/******************************************************************************
rnd()

Pick a random integer between low and high
*****************************************************************************/ 
int NSGAII::rnd(int low, int high) 
{ 
   int i; 
 
   if(low >= high) 
   {
      i = low; 
   }
   else 
   { 
      i = (int)((randomperc() * (high - low + 1)) + low); 
      if(i > high) 
      {
         i = high;
      } 
   } 
   return(i); 
} /* end rnd() */
 
/******************************************************************************
rndreal()

Pick a real random number between specified limits
*****************************************************************************/ 
float NSGAII::rndreal(float lo , float hi) 
{ 
   return((randomperc() * (hi - lo)) + lo); 
} /* end rndreal() */

/******************************************************************************
ranking()

demarcate the different Pareto Fronts
*****************************************************************************/ 
void NSGAII::ranking(NSGAII_Population *pop_ptr)
{
   int i,j,k; /* counters */
   int rnk; /* rank */
   int val; /* value obtained after comparing two individuals */
   int nondom; /* no of non dominated members */
   int maxrank1; /* Max rank of the population */
   int * rankarr; /* [maxpop] Array storing the individual number at a rank */
   int q;
  
   float *ptr1,*ptr2;

   rankarr = (int *)malloc(sizeof(int) * m_maxpop);
  
   /*----------------------------------------------------------------
   RANKING 
   ----------------------------------------------------------------*/
  
   /*Initializing the ranks to zero*/
   rnk = 0 ;
  
   nondom = 0 ;
   maxrank1 = 0;
  
   /*----------------------------------------------------------------
   min_fit is initialize to start distributing the dummy fitness = 
   popsize to the rank one individuals and keeping the record such 
   that the minimum fitness of the better rank individual is always 
   greater than max fitness of the relatively worse rank
   ----------------------------------------------------------------*/
  
   /*----------------------------------------------------------------
   Difference in the fitness of minimum dummy fitness of better rank 
   and max fitness of the next ranked individuals
   ----------------------------------------------------------------*/
  
   /* Initializing all the flags to 2 */
   for(j = 0 ; j < m_popsize; j++)
   {
      pop_ptr->ind[j].flag = 2;
   }
  
   q = 0;
  
   for(k =  0; k < m_popsize; k++,q=0)
   {
      for(j = 0;j  < m_popsize; j++)
      {
         if (pop_ptr->ind[j].flag != 1)
         {
            break; 
         }
      }/* end for() */
      /* Break if all the individuals are assigned a rank */
      if(j == m_popsize)
      {
         break;
      }
      
      rnk = rnk + 1;
      
      for(j = 0 ; j < m_popsize; j++)
      {
         if(pop_ptr->ind[j].flag == 0)
         {
            pop_ptr->ind[j].flag = 2; /* Set the flag of dominated individuals to 2 */
         }
      }/* end for() */
      
      for(i = 0; i < m_popsize; i++)
      {
         /*Select an individual which rank to be assigned*/
	  
         pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
	  
         if(pop_ptr->ind_ptr->flag != 1 && pop_ptr->ind_ptr->flag != 0) 
         {
            ptr1 = &(pop_ptr->ind_ptr->fitness[0]);
	      
            for(j = 0; j < m_popsize ; j++)
            {		  
               /*Select the other individual which has not got a rank*/
               if(i != j)
               {
                  if(pop_ptr->ind[j].flag != 1)
                  {
                     pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
                     ptr2 = &(pop_ptr->ind_ptr->fitness[0]);
			  
                     /*Compare the two individuals for fitness*/
                     val = indcmp(ptr1,ptr2);
			  
                     /*VAL = 2 for dominated individual which rank to be given */
                     /*VAL = 1 for dominating individual which rank to be given */			  
                     /*VAL = 3 for non comparable individuals */			  
                     if( val == 2)
                     { 
                        pop_ptr->ind[i].flag = 0;/* individual 1 is dominated */
                        break;
                     }			  
                     if(val == 1)
                     {
                        pop_ptr->ind[j].flag = 0;/* individual 2 is dominated */
                     }		  
                     if(val == 3)
                     {
                        nondom++;/* individual 1 & 2 are non dominated */
                        if(pop_ptr->ind[j].flag != 0)
                        {
                           pop_ptr->ind[j].flag = 3;
                        }
                     }
                  }   /*end if() */
               } /* end if(i != j) */
            }  /* end for() */
            if(j == m_popsize)
            { 
               /*Assign the rank and set the flag*/
               pop_ptr->ind[i].rank = rnk;
               pop_ptr->ind[i].flag = 1;
               rankarr[q] = i;
               q++;
            }
         }  /* end if(flag check) */
      }/* end for() */
      pop_ptr->rankno[rnk-1] = q;
   }/* end for() */
   maxrank1 = rnk;
  
   /* Find Max Rank of the population */
   for(i = 0; i < m_popsize; i++)
   {
      rnk = pop_ptr->ind[i].rank;
      
      if(rnk > maxrank1) 
      { 
         maxrank1 = rnk;
      }
   }/* end for() */
  
   pop_ptr->maxrank = maxrank1;

   free(rankarr);  
}/* end ranking() */

/******************************************************************************
indcmp()

Compare two individuals
*****************************************************************************/ 
int NSGAII::indcmp(float *ptr1,float *ptr2)
{
   float * fit1; /* [maxfun] */
   float * fit2; /* [maxfun] */
   int i,value,m,n;

   fit1 = (float *)malloc(sizeof(float) * m_maxfun);
   fit2 = (float *)malloc(sizeof(float) * m_maxfun);

   for(i = 0; i < m_nfunc; i++)
   {
      fit1[i] = *ptr1++;
      fit2[i] = *ptr2++;
   }
   m = 0;
   n = 0;
   while(m < m_nfunc && fit1[m] <= fit2[m]) 
   {
      if((fit2[m] -  fit1[m]) < 1e-7) 
      {
         n++;
      }
      m++;
   }/* end while() */
   if(m == m_nfunc) 
   {
      if(n == m_nfunc) 
      {
         value = 3;
      }
      else 
      {
         value = 1; /* value = 1 for dominationg */
      }
   }/* end if() */
   else 
   {
      m = 0;
      n = 0;
      while(m < m_nfunc && fit1[m] >= fit2[m])
      {
	      if((fit1[m] - fit2[m]) < 1e-7) 
         {
            n++;
         }
	   m++;
      }/* end while() */

      if(m == m_nfunc)
      {
	      if(n != m_nfunc)
         {
	         value = 2; /*value =  2 for dominated */
	      }
         else
         {
            value = 3;
         }
      }/* end if() */
      else
      {
         value = 3; /*value = 3 for incomparable */
      }
  }/* end else() */

   free(fit1);
   free(fit2);
  
   return value;
}/* end indcmp() */

/******************************************************************************
realcross()

Crossover for real-coded variables
*****************************************************************************/ 
void NSGAII::realcross(NSGAII_Population *new_pop_ptr,NSGAII_Population *mate_pop_ptr)
{
   int i,j,y,n;
   // int r;
   float rnd,par1,par2,chld1,chld2,betaq,beta,alpha;
   float y1,y2,yu,yl,expp;

   y=0; 
   n=0;
   for(i = 0; i < m_popsize/2; i++)
   {
      rnd = randomperc();
      
      /*Check Whether the cross-over to be performed*/
      if(rnd <= m_pcross)
      {
	      /*Loop over no of variables*/
         for(j = 0; j < m_nvar; j++)
         { 
            /*Selected Two Parents*/ 
            par1 = mate_pop_ptr->ind[y].xreal[j];
            par2 = mate_pop_ptr->ind[y+1].xreal[j]; 
	      
            yl = m_lim_r[j][0];
            yu = m_lim_r[j][1];
	      
            rnd = randomperc();
	      
            /* Check whether variable is selected or not*/
            if(rnd <= 0.5)
            {
               /*Variable selected*/
               m_ncross++;
		  
               if(fabs(par1 - par2) > 0.000001) // changed by Deb (31/10/01)
               {
                  if(par2 > par1)
                  {
                     y2 = par2;
                     y1 = par1;
                  }
                  else
                  {
                     y2 = par1;
                     y1 = par2;
                  }
		      
                  /*Find beta value*/
                  if((y1 - yl) > (yu - y2))
                  {
                     beta = 1 + (2*(yu - y2)/(y2 - y1));
                  }
                  else
                  {
                     beta = 1 + (2*(y1-yl)/(y2-y1));
                  }
		      
                  /*Find alpha*/
                  expp = (float)(m_di + 1.0);
		      
                  beta = (float)(1.0/beta);
		      
                  alpha = (float)(2.0 - pow(beta,expp));
		      
                  if (alpha < 0.0) 
                  {
                     printf("ERRRROR %f %d %d %f %f\n",alpha,y,n,par1,par2);
                     exit(-1);
                  }

                  rnd = randomperc(); 
		      
                  if (rnd <= 1.0/alpha)
                  {
                     alpha = alpha*rnd;
                     expp = (float)(1.0/(m_di + 1.0));
                     betaq = (float)pow(alpha,expp);
                  }/* end if() */
                  else
                  {
                     alpha = alpha*rnd;
                     alpha = (float)(1.0/(2.0-alpha));
                     expp = (float)(1.0/(m_di + 1.0));
                     if (alpha < 0.0) 
                     {
                        printf("ERRRORRR \n");
                        exit(-1);
                     }
                     betaq = (float)pow(alpha,expp);
                  }/* end else() */
		      
                  /*Generating two children*/
                  chld1 = (float)(0.5*((y1+y2) - betaq*(y2-y1)));
                  chld2 = (float)(0.5*((y1+y2) + betaq*(y2-y1)));		      
               }/* end if() */
               else
               {
                  betaq = 1.0;
                  y1 = par1; y2 = par2;
		      
                  /*Generation two children*/
                  chld1 = (float)(0.5*((y1+y2) - betaq*(y2-y1)));
                  chld2 =  (float)(0.5*((y1+y2) + betaq*(y2-y1)));		      
               }/* end else() */
               // added by deb (31/10/01)
               if (chld1 < yl) 
               {
                  chld1 = yl;
               }
               if (chld1 > yu) 
               {
                  chld1 = yu;
               }
               if (chld2 < yl) 
               {
                  chld2 = yl;
               }
               if (chld2 > yu) 
               {
                  chld2 = yu;
               }
            }/* end if() */
            else
            {		  
               /*Copying the children to parents*/
               chld1 = par1;
               chld2 = par2;
            }/* end else() */
            new_pop_ptr->ind[n].xreal[j] = chld1;
            new_pop_ptr->ind[n+1].xreal[j] = chld2;
         }/* end for() */
      }/* end if() */
      else
      {
         for(j = 0; j < m_nvar; j++)
         {
            par1 = mate_pop_ptr->ind[y].xreal[j];
            par2 = mate_pop_ptr->ind[y+1].xreal[j]; 
            chld1 = par1;
            chld2 = par2;
            new_pop_ptr->ind[n].xreal[j] = chld1;
            new_pop_ptr->ind[n+1].xreal[j] = chld2;
         }/* end for() */
      }/* end else() */
      n = n+2; 
      y=y+2;
   }/* end for() */
}/* end realcross() */

/******************************************************************************
realinit()

Initialize the real-coded population
*****************************************************************************/ 
void NSGAII::realinit(NSGAII_Population *pop_ptr)
{
   int i,j;
   float d,d1;
  
   for (i = 0 ; i < m_popsize ; i++)
   { 
      for (j = 0; j < m_nvar; j++)
      {
         d = randomperc();
         d1 = 2*d - 1;
         /*if limits are not specified then generate any number between zero and infinity*/
         if(m_ans != 1)
         {
            pop_ptr->ind[i].xreal[j] = 1/d1;
         }
	  
         /*if limits are specified it generates the value in range of minimum and maximum value of the variable*/
         else
         {
            pop_ptr->ind[i].xreal[j] = d*(m_lim_r[j][1] - m_lim_r[j][0])+m_lim_r[j][0];
         }
      }/* end for() */
   }/* end for() */
}/* end realinit() */

/******************************************************************************
real_mutate()
*****************************************************************************/ 
void NSGAII::real_mutate(NSGAII_Population *new_pop_ptr)
{
   int i,j;
   float rnd,delta,indi,deltaq;
   float y,yl,yu,val,xy;
  
   for(j = 0; j < m_popsize; j++)
   {
      for (i = 0; i < m_nvar; i++)
      {
         rnd = randomperc();
	  
         /*For each variable find whether to do mutation or not*/
         if(rnd <= m_pmut_r)
         {
            y = new_pop_ptr->ind[j].xreal[i];
            yl = m_lim_r[i][0];
            yu = m_lim_r[i][1];
	      
            if(y > yl)
            {
               /*Calculate delta*/
		  
               if((y-yl) < (yu-y))
               {
                  delta = (y - yl)/(yu - yl);
               }
               else
               {
                  delta = (yu - y)/(yu-yl);
		         }

               rnd = randomperc(); 
		  
               indi = (float)(1.0/(m_dim +1.0));
		  
               if(rnd <= 0.5)
               {
                  xy = (float)(1.0-delta);
                  val = (float)(2*rnd+(1-2*rnd)*(pow(xy,(m_dim+1))));
                  deltaq =  (float)(pow(val,indi) - 1.0);
               }
               else
               {
                  xy = (float)(1.0-delta);
                  val = (float)(2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(m_dim+1))));
                  deltaq = (float)(1.0 - (pow(val,indi)));
               }
		  
               /*Change the value for the parent */
               //  *ptr  = *ptr + deltaq*(yu-yl);
               // Added by Deb (31/10/01)
               y = y + deltaq * (yu-yl);
               if (y < yl) 
               { 
                  y=yl;
               } 
               if (y > yu) 
               { 
                  y=yu;
               }
               new_pop_ptr->ind[j].xreal[i] = y;
            }/* end if() */
            else // y == yl 
            {
               xy = randomperc();
               new_pop_ptr->ind[j].xreal[i] = xy*(yu - yl) + yl;
            }
            m_nmut++;
         }/* end if() */
      }/* end for() */
   }/* end for() */
}/* end real_mutate() */

/******************************************************************************
rselect()

selection operator for real-coded variables
*****************************************************************************/ 
void NSGAII::rselect(NSGAII_Population *old_pop_ptr,NSGAII_Population *mate_pop_ptr)
{
   int *fit_ptr1,*fit_ptr2;
  
   float rnd2,*f1_ptr,*f2_ptr;
  
   float *s1_ptr,*s2_ptr;
  
   float *select_ptr;
  
   NSGAII_Individual *j,*j1;
  
   int i,rnd,rnd1,n,j2,r,s;
  
   old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
   mate_pop_ptr->ind_ptr= &(mate_pop_ptr->ind[0]); 
   j =  &(old_pop_ptr->ind[m_popsize - 1]);
  
   old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]); 
   j2 = 0;
   r = m_popsize;
   s = m_chrom;

   for(n = 0; n < m_popsize; n++)
   {
      mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[n]);
      select_ptr = &(mate_pop_ptr->ind_ptr->gener[0]);
      
      rnd2 = randomperc(); 
      
      rnd2 = m_popsize * rnd2; 
      
      rnd = (int)floor(rnd2);
      if(rnd == 0)
      rnd = m_popsize - n; // -2;
      
      if(rnd >=  m_popsize)
      {
         rnd  = (m_popsize - 2)/2; // 4)/2;
      }

      /*Select first individual randomly*/	
      j = &(old_pop_ptr->ind[rnd-1]);
      
      rnd2 = randomperc(); 
      rnd2 = m_popsize * rnd2; 
      rnd1 = (int)floor(rnd2);
      
      if (rnd1 == 0)
      {
         rnd1 = m_popsize - n;
      }

      if(rnd1 >= m_popsize)
      {
         rnd1 = (m_popsize - 4)/2; // 2)/2;
      }

      /*Select second parent randomly */
      j1 = &(old_pop_ptr->ind[rnd1-1]); 
      
      old_pop_ptr->ind_ptr = j; 
      
      /*Assigning Pointers to the genes and dummyfitness */
      s1_ptr = &(old_pop_ptr->ind_ptr->gener[0]); 
      fit_ptr1 = &(old_pop_ptr->ind_ptr->rank); 
      f1_ptr = &(old_pop_ptr->ind_ptr->cub_len); 
       
      old_pop_ptr->ind_ptr = j1; 
      s2_ptr = &(old_pop_ptr->ind_ptr->gener[0]); 
      fit_ptr2 = &(old_pop_ptr->ind_ptr->rank); 
      f2_ptr = &(old_pop_ptr->ind_ptr->cub_len); 

      /*--------------------------------------------------------------------------
      SELECTION PROCEDURE
      --------------------------------------------------------------------------*/
      
      /*Selecting one parent on the basis of tournament selection*/
      if(*fit_ptr1 > *fit_ptr2)
      {
         for(i = 0; i < m_chrom; i++)
         {
            *select_ptr++ = *s2_ptr++;
         }
      }/* end if() */
      else
      {
         if(*fit_ptr1 < *fit_ptr2)
         {
            for(i = 0; i < m_chrom; i++)
            {
               *select_ptr++ = *s1_ptr++;
            }
         }
         else
         {
            if(*f1_ptr < *f2_ptr)
            {
               for(i = 0; i < m_chrom; i++)
               {
                  *select_ptr++ = *s2_ptr++;
               }
            }
            else
            {
               for(i = 0; i < m_chrom; i++)
               {
                  *select_ptr++ = *s1_ptr++;
               }
            }/* end else */
         }/* end else() */
      }/* end else() */
   }/* end for() */
  
   for(i = 0; i < m_popsize; i++)
   {
      // mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[i]);
      for(j2 = 0; j2 < m_chrom; j2++)
      {
         mate_pop_ptr->ind[i].xreal[j2] = mate_pop_ptr->ind[i].gener[j2];
      }/* end for() */
   }/* end for() */
}/* end rselect() */

/******************************************************************************
report()

print out results to file
*****************************************************************************/ 
void NSGAII::report(int t,NSGAII_Population *pop1_ptr,NSGAII_Population *pop2_ptr,FILE *rep_ptr,FILE *gen_ptr, FILE *lastit )
{
   int i,j,*rptr,*rptr1; 

   float *ptr1,*fptr,*fptr1, *ptr1_b, *ptr2_b;

   float *ptr2,*cons_ptr1,*cons_ptr2, *err2;
  
   fprintf(rep_ptr,"\n\n---------------------------------------------------\n");
   fprintf(rep_ptr,"Generation No.     ->%d\n",t+1);
   fprintf(rep_ptr,"------------------------------------------------------\n");
   if(m_ncons == 0)
   {
      fprintf(rep_ptr," variables (real %d binary %d)  fitness (%d)  rank cublen || variables  fitness rank cublen\n",m_nvar, m_nchrom, m_nfunc);
   }
   else
   {
      fprintf(rep_ptr," variables (real %d binary %d)  fitness (%d) constraint (%d) penalty rank cublen || variables  fitness constraint penalty rank cublen\n",m_nvar, m_nchrom, m_nfunc, m_ncons);
   }

   pop1_ptr->ind_ptr = &(pop1_ptr->ind[0]);
  
   pop2_ptr->ind_ptr = &(pop2_ptr->ind[0]); // Deb 31/10/01


   for(i = 0; i < m_popsize; i++)
   {
      fprintf(rep_ptr,"\n------------------------------------------------\n"); 

      ptr1_b = &(pop1_ptr->ind_ptr->xbin[0]);
      ptr2_b = &(pop2_ptr->ind_ptr->xbin[0]);

      ptr1 = &(pop1_ptr->ind_ptr->xreal[0]);
      ptr2 = &(pop2_ptr->ind_ptr->xreal[0]);  // Deb 31/10/01
      
      fptr = &(pop1_ptr->ind_ptr->fitness[0]);
      fptr1 = &(pop2_ptr->ind_ptr->fitness[0]);
      
      rptr = &(pop1_ptr->ind_ptr->rank);
      rptr1 = &(pop2_ptr->ind_ptr->rank);
      
      cons_ptr1 = &(pop1_ptr->ind_ptr->constr[0]);
      cons_ptr2 = &(pop2_ptr->ind_ptr->constr[0]);
      
      err2 = &(pop2_ptr->ind_ptr->error);

      for(j = 0; j < m_nvar; j++)
      {
         fprintf(rep_ptr,"%f ",*ptr1++);
      }

      for(j = 0; j < m_nchrom; j++)
      {
         fprintf(rep_ptr,"%f ",*ptr1_b++);
      }
      if (t == (m_gener - 1))
      {
         for(j = 0; j < m_nfunc; j++)
         {   
            if ((*err2 <= 0.0) && (*rptr1 == 1))
            {
               fprintf(lastit,"%f\t",*fptr1++);
            }
            else 
            { 
               *fptr1++;
            }
         }/* end for() */
         if ((*err2 <= 0.0) && (*rptr1 == 1))
         {
            fprintf(lastit,"\n");
         }
      }/* end if() */
      fptr =  &(pop1_ptr->ind_ptr->fitness[0]);
      fptr1 = &(pop2_ptr->ind_ptr->fitness[0]);
      
      for(j = 0; j < m_nfunc; j++)
      {
         fprintf(rep_ptr,"  %.4f",*fptr++);
      }
      if(m_ncons != 0)
      {
         for(j = 0; j < m_ncons;j++)
         {
            fprintf(rep_ptr,"  %.2e",*cons_ptr1++);
         }
         fprintf(rep_ptr," %.2e",pop1_ptr->ind_ptr->error);
      }/* end if() */
      
      fprintf(rep_ptr," %d ",*rptr);
      
      fprintf(rep_ptr,"%f ",pop1_ptr->ind_ptr->cub_len);
      fprintf(rep_ptr,"|**|");

      for(j = 0; j < m_nvar; j++)
      {
         fprintf(rep_ptr," %f ",*ptr2);
         fprintf(gen_ptr," %f ",*ptr2++);
      }
      for(j = 0; j < m_nchrom; j++)
      {
         fprintf(rep_ptr,"%f ",*ptr2_b); 
         fprintf(gen_ptr,"%f ",*ptr2_b++);
      }
      for(j = 0; j < m_nfunc; j++)
      {	
         fprintf(rep_ptr,"  %f",*fptr1);
         fprintf(gen_ptr,"  %f",*fptr1++);
      }
      fprintf(rep_ptr," %d ",*rptr1);
      
      if(m_ncons != 0)
      {
         for(j = 0; j < m_ncons; j++)
         {
            fprintf(rep_ptr,"  %.2e",*cons_ptr2);
            fprintf(gen_ptr,"  %.2e",*cons_ptr2++);
         }/* end for() */
         fprintf(rep_ptr," %.2e",pop2_ptr->ind_ptr->error);
         fprintf(gen_ptr," %.2e",pop2_ptr->ind_ptr->error);
      }/* end if() */
      
      fprintf(rep_ptr," %f ",pop2_ptr->ind_ptr->cub_len);
      
      pop1_ptr->ind_ptr = &(pop1_ptr->ind[i+1]);
      pop2_ptr->ind_ptr = &(pop2_ptr->ind[i+1]);
      
      fprintf(gen_ptr,"\n");
   }/* end for() */
  
   fprintf(rep_ptr,"\n--------------------------------------------------\n\n"); 
   fprintf(rep_ptr,"-------------------------------------------------------\n");
   fprintf(gen_ptr,"\n--------------------------------------------------\n\n");
}/* end report() */

/******************************************************************************
roulette()

roulette wheel selection operator
*****************************************************************************/ 
void NSGAII::roulette(NSGAII_Population *old_pop_ptr,NSGAII_Population *mate_pop_ptr)
{
   int i,j,k,r,*gene1,*gene2;
   float sum,c,rnd,*fit;
   float * cumprob; /* [maxpop] */ 

   cumprob = (float *)malloc(sizeof(float)*m_maxpop);

   /*Find Sum of the fitness*/
   sum = 0;
  
   for (i = 0; i < m_popsize; i++)
   { 
      old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[i]);
      fit = &(old_pop_ptr->ind_ptr->dumfit);
      r = old_pop_ptr->ind_ptr->rank;
      sum = sum +  (*fit) * 1000 * (m_popsize - r);
   }/* end for() */

   c = 0;

   /*Finding the cumulative probability for the different individuals*/
   for(i = 0; i < m_popsize; i++)
   {
      old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[i]);
      
      fit = &(old_pop_ptr->ind_ptr->dumfit);
      
      r = old_pop_ptr->ind_ptr->rank;
      
      c = c +  (*fit)* 1000 * (m_popsize - r);
      
      cumprob[i] = c/sum ;
   }/* end for() */
  
   for(i = 0; i < m_popsize; i++)
   {
      rnd = randomperc();
      j = 0;
     
      /*Generate a Random No to find which individual to select*/  
      while(rnd > cumprob[j])
      {
         j++;
      }
     
      old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
      mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[i]);
     
      if(m_optype != 3)
      {
         gene1 = &(old_pop_ptr->ind_ptr->genes[0]);
	 
         gene2 = &(mate_pop_ptr->ind_ptr->genes[0]);
	 
         for(k = 0; k < m_chrom; k++)
         {
            *gene2++ = *gene1++;
         }	 
      }/* end if() */
      else
      {
         for(k = 0; k < m_chrom; k++)
         {
            mate_pop_ptr->ind_ptr->gener[k] = old_pop_ptr->ind_ptr->gener[k];
         }
      }/* end else() */
      for(k = 0; k < m_nvar; k++)
      {
         mate_pop_ptr->ind_ptr->xreal[k] = old_pop_ptr->ind_ptr->xreal[k];
      }
   }/* end for() */
   free(cumprob);
}/* end roulette() */

/******************************************************************************
nselect()

selection operator for binary coded variables
*****************************************************************************/ 
void NSGAII::nselect(NSGAII_Population *old_pop_ptr,NSGAII_Population *pop2_ptr)
{
   int *fit_ptr1,*fit_ptr2;

   float rnd2,*f1_ptr,*f2_ptr;
  
   int *s1_ptr,*s2_ptr,*select_ptr;

   float *select_ptr_r, *s1_ptr_r, *s2_ptr_r;
  
   NSGAII_Individual *j,*j1;
  
   int i,rnd,rnd1,k,n,j2,r,s;
  
   old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
  
   pop2_ptr->ind_ptr= &(pop2_ptr->ind[0]); 
  
   j =  &(old_pop_ptr->ind[m_popsize-1]);
  
   old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]); 
   j2 = 0;
   r = m_popsize;
   s = m_chrom;
  
   for(n = 0,k = 0; n < m_popsize; n++,k++)
   {
      pop2_ptr->ind_ptr = &(pop2_ptr->ind[k]);
      select_ptr = &(pop2_ptr->ind_ptr->genes[0]);
      select_ptr_r = &(pop2_ptr->ind_ptr->xreal[0]);

      rnd2 = randomperc(); 

      rnd2 = m_popsize * rnd2; 

      rnd = (int)floor(rnd2);

      if(rnd == 0)
      {
         rnd = m_popsize - k;
      }
      if(rnd == m_popsize)
      {
         rnd = (m_popsize-2)/2;
      }
      /*Select first parent randomly*/	
      j = &(old_pop_ptr->ind[rnd-1]);
      
      rnd2 = randomperc(); 
      
      rnd2 = m_popsize * rnd2; 
      
      rnd1 = (int)floor(rnd2);
      
      if (rnd1 == 0)
      {
         rnd1 = m_popsize - n;
      }

      if(rnd1 == m_popsize)
      {
         rnd1 = (m_popsize - 4)/2;
      }
      
      /*Select second parent randomly*/
      j1 = &(old_pop_ptr->ind[rnd1-1]);
      
      old_pop_ptr->ind_ptr = j;
      
      s1_ptr = &(old_pop_ptr->ind_ptr->genes[0]);
      s1_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
      fit_ptr1 = &(old_pop_ptr->ind_ptr->rank);
      f1_ptr = &(old_pop_ptr->ind_ptr->cub_len);
      
      old_pop_ptr->ind_ptr = j1;
      s2_ptr = &(old_pop_ptr->ind_ptr->genes[0]);
      s2_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
      fit_ptr2 = &(old_pop_ptr->ind_ptr->rank);
      f2_ptr = &(old_pop_ptr->ind_ptr->cub_len);

      /*--------------------------------------------------------------------------
      SELECTION PROCEDURE
      --------------------------------------------------------------------------*/
      
      /*Comparing the fitnesses*/
      if(*fit_ptr1 > *fit_ptr2)
      {
         for(i = 0; i < m_chrom; i++)
         {
            *select_ptr++=*s2_ptr++;
         }
         for(i = 0; i < m_nvar; i++)
         {
            *select_ptr_r++=*s2_ptr_r++;
         }
      }/* end if() */
      else
      {
         if(*fit_ptr1 < *fit_ptr2)
         {
            for(i = 0; i < m_chrom; i++)
            {
               *select_ptr++ = *s1_ptr++;
            }
            for(i = 0; i < m_nvar; i++)
            {
               *select_ptr_r++ = *s1_ptr_r++;
            }
         }/* end if() */
         else
         {
            if(*f1_ptr < *f2_ptr)
            {
               for(i = 0; i < m_chrom; i++)
               {
                  *select_ptr++ = *s2_ptr++;
               }
               for(i = 0;i < m_nvar; i++)
               {
                  *select_ptr_r++ = *s2_ptr_r++;
               }
            }/* end if() */
            else
            {
               for(i = 0; i < m_chrom; i++)
               {
                  *select_ptr++ = *s1_ptr++;
               }
               for(i = 0; i < m_nvar; i++)
               {
                  *select_ptr_r++ = *s1_ptr_r++;
               }
            }/* end else() */
         }/* end else() */
      }/* end else() */
   }/* end for() */
}/* end nselect() */

/******************************************************************************
unicross()

uniform crossover operator
*****************************************************************************/ 
void NSGAII::unicross(NSGAII_Population *new_pop_ptr, NSGAII_Population *mate_pop_ptr)
{
   int i,j,*gene,y,n,*par1,*par2,*chld1,*chld2;
   float rnd;

   for(i = 0,y = 0,n = 0; i < m_popsize; i++)
   {
      for(j = 0; j < m_chrom; j++)
      {
         /*Select a bit for doing cross-over*/	
         new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[y]);
         chld1 = &(new_pop_ptr->ind_ptr->genes[j]);

         new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[y+1]);
         chld2 = &(new_pop_ptr->ind_ptr->genes[j]);
	
         mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[n]);
         par1 = &(mate_pop_ptr->ind_ptr->genes[j]);

         mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[n+1]);
         par2 = &(mate_pop_ptr->ind_ptr->genes[j]);

         rnd = randomperc();

         /* Checking whether to do cross-over or not */
         if(rnd <= m_pcross)
         {
            m_ncross++;
            *chld1 = *par2;
            *chld2 = *par2;
         }
         else
         {
            *chld1 = *par1;
            *chld2 = *par2;
         }
      }/*end for() */

      y = y+2;
      n = n+2;
   }/* end for() */

   for(i = 0; i < m_popsize; i++)
   {
      new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[i]);
      gene = &(new_pop_ptr->ind_ptr->genes[0]);
      for(j = 0; j < m_chrom; j++)
      {
         gene = &(new_pop_ptr->ind_ptr->genes[j]);
      }/* end for() */
   }/* end for() */
}/* end unicross() */

/******************************************************************************
NSGAII_Program()

Calibrate or optimize the model using NSGAII.
******************************************************************************/
void NSGAII_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("NSGAII", 1);
   NSGAII * TheAlg = new NSGAII(model);
   MEM_CHECK(TheAlg);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end NSGAII_Program() */




