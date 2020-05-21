/******************************************************************************
File      : SCEUA.cpp
Author    : L. Shawn Matott - Converted from Origincal SCEUA Fortran code
Copyright : 2009, L. Shawn Matott

The SCE-UA method is a general purpose global optimization
program - Shuffled Complex Evolution.  It was originally developed 
by Dr. Qingyun Duan as part of his doctoral dissertation work at
University of Arizona, Tucson, AZ 85721, USA. 

The dissertation is entitled "A Global Optimization Strategy for
Efficient and Effective Calibration of Hydrologic Models".  The
program has since been modified to make it easier for use on
problems of users' interests.  

The algorithm has been described in detail in an article entitled :
"Effective and Efficient Global Optimization for Conceptual Rainfall-Runoff Models", 
Water Resources Research, Vol 28(4), pp.1015-1031, 1992

And in an article entitled:
"A Shuffled Complex Evolution Approach for Effective and Efficient Global Minimization",
 by Q. Duan, V.K. Gupta and S. Sorooshian, Journal of Optimization Theory and its 
Applications, Vol 76(3), pp 501-521, 1993.  

A paper entitled "Optimal Use of the SCE-UA Global Optimization Method for Calibrating Watershed Models", 
by Q. Duan, S. Sorooshian, & V.K. Gupta, Journal of Hydrology, Vol.158, 265-284, 1994, 
discussed how to use the SCE-UA Method in an efficient and  effective manner.

Input Summary for the SCEUA algorithm (adapted from original Fortran-based
description):
==========================================================================
variable   type     description
MAXN       integer  Maximum number of trials allowed before
                    optimization is terminated.  The purpose of
                    MAXN is to stop an optimization search before
                    too much computer time is expended.  MAXN
                    should be set large enough so that optimization
                    is generally completed before MAXN trials are
                    performed. Recommended value is 10,000 (increase or
                    decrease as necessary).
---> this parameter is called m_Budget within Ostrich

KSTOP      integer  Number of shuffling loops in which the 
                    criterion must improve by the specified
                    percentage or else optimization will be
                    terminated. Recommended value: 5.
---> this parameter is called m_Kstop within Ostrich

PCENTO     double   Percentage by which the criterion value must
                    change in the specified number of shuffling 
                    loops or else optimization is terminated
                    (Use decimal equivalent: Percentage/100).
                    Recommended value: 0.01.
---> this parameter is called m_Pcento within Ostrich

NGS        integer  Number of complexes used for optimization
                    search.  Minimum value is 1.
                    Recommended value is between 2 and 20 depending
                    on the number of parameters to be optimized and
                    on the degree of difficulty of the problem.
---> this parameter is called m_NumComplexes within Ostrich

ISEED      integer  Random seed used in optimization search.  Enter
                    any integer number.  Default value (=1969) is
                    assumed if this field is left blank.
                    Recommended value: any large integer.
---> this parameter is called m_Seed within Ostrich

IDEFLT     integer  Flag for setting the control variables of the
                    SCE-UA algorithm.  Enter false or leave the field
                    blank for default setting.  Enter true for user
                    specified setting.
                    If option true is chosen, user must specify alg. 
                    parameters.
---> this parameter is called m_bUserConfig within Ostrich

NPG        integer  Number of points in each complex.  NPG should
                    be greater than or equal to 2.  The default
                    value is equal to (2 * number of optimized
                    parameters + 1).
---> this parameter is called m_PtsPerComplex within Ostrich

NPS        integer  Number of points in each sub-complex.  NPS
                    should be greater than or equal to 2 and less
                    than NPG.  The default value is equal to 
                    (number of optimized parameters + 1).
---> this parameter is called m_PtsPerSubComplex within Ostrich

NSPL       integer  Number of evolution steps taken by each complex
                    before next shuffling.  Default value is equal
                    to NPG.
---> this parameter is called m_NumEvoSteps within Ostrich

MINGS      integer  Minimum number of complexes required for
                    optimization search, if the number of complexes
                    is allowed to reduce as the optimization search
                    proceeds.  The default value is equal to NGS.
---> this parameter is called m_MinComplexes within Ostrich

INIFLG     integer  Flag on whether to include an initial point in
                    the starting population.  Enter true if the initial 
                    point is to be included.  The default value is equal to false.
---> this parameter is called m_bUseInitPt within Ostrich

IPRINT    integer   Print-out control flag.  Enter '0' to print out
                    the best estimate of the global optimum at the
                    end of each shuffling loop.  Enter '1' to print
                    out every point in the entire sample population
                    at the end of each shuffling loop.  The default
                    value is equal to 0. Enter 2 to ignore this variable
                    and use conventional Ostrich output.
---> this parameter is called m_OutputMode within Ostrich

PARAMS     double * Initial estimates of the parameters to be optimized.
---> this parameter is called m_pParams within Ostrich

LOWER      double * Lower bounds of the parameters to be optimized.
---> this parameter is called m_pLower within Ostrich

UPPER      double * Upper bounds of the parameters to be optimized.
---> this parameter is called m_pUpper within Ostrich

Version History
10-31-09    lsm   Created
******************************************************************************/
#include <math.h>
#include <string.h>

#include "SCEUA.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "StatsClass.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
void SCEUA::WarmStart(void)
{
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pModel->GetParamGroupPtr()->WriteParams(pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
}/* end WarmStart() */

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
SCEUA::SCEUA(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pParams = NULL;
   m_pUpper = NULL;
   m_pLower = NULL;
   m_bUseInitPt = false;
   m_fSaved = NEARLY_HUGE;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()
******************************************************************************/
void SCEUA::Destroy(void)
{ 
   delete [] m_pParams;
   delete [] m_pLower;
   delete [] m_pUpper;
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using SCEUA.
******************************************************************************/
void SCEUA::Calibrate(void)
{ 
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);
   
   Optimize();

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   id = 0;
   sprintf(fileName, "OstOutput%d.txt", id);

   //write statistics of best parameter set to output file
   pFile = fopen(fileName, "a");   
   m_pStats->WriteStats(pFile);
   fclose(pFile);

   //write statistics of best parameter set to output file
   m_pStats->WriteStats(stdout);
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using SCEUA.
******************************************************************************/
void SCEUA::Optimize(void)
{
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();

   InitFromFile(GetInFileName());

   WriteSetup(m_pModel, "Shuffled Complex Evolution - University of Arizona");
   m_CurIter = 0;
   //write banner
   WriteBanner(m_pModel, "gen   best value     ", "Pct. Complete");
  
   scemain(); //main SCE implemenation, converted from FORTRAN

   //place model at optimal prameter set
   pGroup->WriteParams(m_pParams);
   m_pModel->Execute();

   WriteOptimal(m_pModel, m_Best);
   m_pStatus.pct = 100.00;
   m_pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&m_pStatus);
   //write algorithm metrics
   WriteAlgMetrics(this);
} /* end Optimize() */

/******************************************************************************
scemain()

THIS IS THE MAIN PROGRAM CALLING SUBROUTINES SCEIN AND SCEUA
******************************************************************************/
void SCEUA::scemain(void)
{
   double * a, * bl, * bu;
   int jseed[10] = {2,3,5,7,11,13,17,19,23,29}; 

   a = m_pParams;
   bl = m_pLower;
   bu = m_pUpper;

   if(m_OutputMode != 2)
      printf(" ENTER THE MAIN PROGRAM --- \n"); 

   int nopt, kstop, iseed, ngs, npg, nps, nspl, mings, iprint, maxn;
   int iniflg;
   double pcento;
   nopt = m_np = m_pModel->GetParamGroupPtr()->GetNumParams();
   scein(a,bl,bu,nopt,&maxn,&kstop,&pcento,&iseed,
         &ngs,&npg,&nps,&nspl,&mings,&iniflg,&iprint);

   int nrun = 1;

   int i;
   for(i = 0; i < nrun; i++)
   {
      if (nrun != 1) iseed = jseed[i];
      if(m_OutputMode != 2)
         printf("@ SCE-UA Run Number %d Random Seed Value %d\n",i,iseed);
      sceua(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
            ngs,npg,nps,nspl,mings,iniflg,iprint);
   }/* end for() */
}/* end scemain() */

/******************************************************************************
scein()

THIS SUBROUTINE READS AND PRINTS THE INPUT VARIABLES FOR
SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
 -- Version 2.1

ADAPTED FROM FORTRAN CODE WRITTEN BY 
    QINGYUN DUAN - UNIVERSITY OF ARIZONA, APRIL 1992
******************************************************************************/
void SCEUA::scein
(
   double * a,
   double * bl,
   double * bu,
   int nopt,
   int *maxn,
   int *kstop,
   double *pcento,
   int *iseed,
   int *ngs,
   int *npg,
   int *nps,
   int *nspl,
   int *mings,
   int *iniflg,
   int *iprint
)
{
   char pcntrl[100],deflt[100],usrsp[100];
   char reduc[40],initl[40],ysflg[40],noflg[40],**xname; 

   strcpy(deflt, "DEFAULT"); 
   strcpy(usrsp, "USER SPEC."); 
   strcpy(ysflg, "YES"); 
   strcpy(noflg, "NO");  
   xname = new char*[nopt];
   for(int i = 0; i < nopt; i++)
   {
     xname[i] = new char[50];   
     strncpy(xname[i], m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetName(), 9);
   }

   if(m_OutputMode != 2) printf("ENTER THE SCEIN SUBROUTINE --- \n");

   //INITIALIZE I/O VARIABLES
   FILE * pIn  = fopen("sce.in", "r");
   FILE * pOut = fopen("sce.out", "w");

   int ierror = 0;
   int iwarn = 0;
   if(m_OutputMode != 2)
   {
      fprintf(pOut, "\
          SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION\n\
          ==============================================\n\n\n");
   }

   // READ THE SCE CONTROL PARAMETERS
   int ideflt = 0;
   char line[160];
   fgets(line, 1600, pIn);
   sscanf(line, "%d %d %lf %d %d %d", maxn, kstop, pcento, ngs, iseed, &ideflt);
   if (*iseed == 0) *iseed = 1969;

  //IF ideflt IS EQUAL TO 1, READ THE SCE CONTROL PARAMETERS
   if (ideflt == 1)
   {
      fgets(line, 1600, pIn);
      sscanf(line, "%d %d %d %d %d %d", npg, nps, nspl, mings, iniflg, iprint);
      strcpy(pcntrl, usrsp);
   }
   else
   {
      strcpy(pcntrl, deflt);
   }

   //READ THE INITIAL PARAMETER VALUES AND THE PARAMETER BOUNDS
   int iopt;
   for(iopt = 0; iopt < nopt; iopt++)
   {
      fgets(line, 1600, pIn);
      sscanf(line,"%lf %lf %lf", &(a[iopt]), &(bl[iopt]), &(bu[iopt]));
   }

   //IF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO THE DEFAULT VALUES
   if (ideflt == 0)
   {
      *npg = 2*nopt + 1;
      *nps = nopt + 1;
      *nspl = *npg;
      *mings = *ngs;
      *iniflg = 0;
      *iprint = 0;
   }/* end if() */

   //CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
   if ((*ngs < 1) || (*ngs >= 1320))
   {
      fprintf(pOut, 
              "**ERROR** NUMBER OF COMPLEXES IN INITIAL POPULATION (%d) IS NOT A VALID CHOICE\n",
              *ngs);
      ierror = ierror + 1;
   }

   if ((*kstop < 0) || (*kstop >= 20))
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF SHUFFLING LOOPS IN \
WHICH THE CRITERION VALUE MUST CHANGE SHOULD BE \
GREATER THAN 0 AND LESS THAN 10. kstop = %d WAS SPECIFIED. \
BUT kstop = 5 WILL BE USED INSTEAD.\n", 
      *kstop);
      iwarn = iwarn + 1;
      *kstop = 5;
   }/* end if() */

   if((*mings < 1) || (*mings > *ngs))
   {
      fprintf(pOut, 
"**WARNING** THE MINIMUM NUMBER OF COMPLEXES (%d) \
IS NOT A VALID CHOICE. SET IT TO DEFAULT \n", 
      *mings);
      iwarn = iwarn + 1;
      *mings = *ngs;
   }/* end if() */

   if ((*npg < 2) || (*npg > 1320/MyMax(*ngs,1))) 
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF POINTS IN A COMPLEX (%d) \
IS NOT A VALID CHOICE, SET IT TO DEFAULT\n", 
      *npg);
      iwarn = iwarn + 1;
      *npg = 2*nopt+1;
   }/* end if() */

   if ((*nps < 2) || (*nps > *npg) || (*nps > 50))
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF POINTS IN A SUB-COMPLEX (%d) \
IS NOT A VALID CHOICE, SET IT TO DEFAULT\n",
      *nps);
      iwarn = iwarn + 1;
      *nps = nopt + 1;
   }/* end if() */

   if (*nspl < 1)
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF EVOLUTION STEPS \
TAKEN IN EACH COMPLEX BEFORE SHUFFLING (%d) \
IS NOT A VALID CHOICE, SET IT TO DEFAULT\n", 
      *nspl);
      iwarn = iwarn + 1;
      *nspl = *npg;
   }/* end if() */
   
   // COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPULATION
   int npt;
   npt = (*ngs) * (*npg); //npt = ngs * npg

   if (npt > 1320)
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF POINTS IN INITIAL \
POPULATION (%d) EXCEED THE POPULATION LIMIT \
SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS\n", 
      npt);
      iwarn = iwarn + 1;
      *ngs = 2;
      *npg = 2*nopt + 1;
      *nps = nopt + 1;
      *nspl = *npg;
   } /* end if() */

   // PRINT OUT THE TOTAL NUMBER OF ERROR AND WARNING MESSAGES
   if (ierror >= 1)
      fprintf(pOut, "*** TOTAL NUMBER OF ERROR MESSAGES IS %d\n",ierror);

   if (iwarn >= 1)
      fprintf(pOut, "*** TOTAL NUMBER OF WARNING MESSAGES IS %d\n",iwarn);

   if (*mings < *ngs)
      strcpy(reduc, ysflg);
   else
      strcpy(reduc, noflg);

   if (iniflg != 0)
      strcpy(initl, ysflg);
   else
      strcpy(initl, noflg);

   //PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
   fprintf(pOut,"\
  SCE CONTROL     MAX TRIALS     REQUIRED IMPROVEMENT     RANDOM\n\
   PARAMETER        ALLOWED      PERCENT    NO. LOOPS      SEED\n\
  -----------     ----------     -------    ---------     ------\n");

   double pcenta;
   pcenta=(*pcento)*100.;
   fprintf(pOut,"  %-11s     %-10d     %7.2lf    %-9d     %-6d\n\n\n",
   pcntrl, *maxn, pcenta, *kstop, *iseed);

   fprintf(pOut,"\
                  SCE ALGORITHM CONTROL PARAMETERS\n\
                  ================================\n\n\
  NUMBER OF     POINTS PER     POINTS IN      POINTS PER    EVOL. STEPS\n\
  COMPLEXES      COMPLEX      INI. POPUL.     SUB-COMPLX    PER COMPLEX\n\
  ---------     ----------    -----------     ----------    -----------\n");
   fprintf(pOut,"  %-9d     %-10d    %-11d     %-10d    %-11d\n\n\n", 
           *ngs, *npg, npt, *nps, *nspl);

   fprintf(pOut,"\
               COMPLX NO.     MIN COMPLEX     INI. POINT\n\
               REDUCTION      NO. ALLOWED      INCLUDED\n\
               ----------     -----------     ----------\n");
   fprintf(pOut,"               %-10s     %-11d     %-10s\n\n\n", 
           reduc, *mings, initl);

   fprintf(pOut,"\
        INITIAL PARAMETER VALUES AND PARAMETER BOUNDS\n\
        =============================================\n\n\
  PARAMETER     INITIAL VALUE     LOWER BOUND     UPPER BOUND\n\
  ---------     -------------     -----------     -----------\n");
   for(int i = 0; i < nopt; i++)
   {
      fprintf(pOut, "  %-9s     %13.3lf     %11.3lf     %11.3lf\n",
              xname[i], a[i], bl[i], bu[i]);
   }/* end for() */
   fprintf(pOut,"\n\n");

   for(int i = 0; i < nopt; i++)
   {
     delete [] xname[i];   
   }
   delete [] xname;

   if (ierror >= 1)
   {
      fprintf(pOut, 
"*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE OF INPUT DATA ERROR ***\n");
      fclose(pIn);
      fclose(pOut);
      ExitProgram(1);
   }/* end if() */

   fclose(pIn);
   fclose(pOut);
}/* end scein() */

/******************************************************************************
sceua()

  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
     -- Version 2.1

  by QINGYUN DUAN (adpated to C++ by L. Shawn Matott)
  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
  (602) 621-9360, email: duan@hwr.arizona.edu

  WRITTEN IN OCTOBER 1990.
  REVISED IN AUGUST 1991
  REVISED IN APRIL 1992

  STATEMENT BY AUTHOR:
  --------------------

     This general purpose global optimization program is developed at
     the Department of Hydrology & Water Resources of the University
     of Arizona.  Further information regarding the SCE-UA method can
     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
     at the address and phone number listed above.  We request all
     users of this program make proper reference to the paper entitled
     'Effective and Efficient Global Optimization for Conceptual
     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.

  LIST OF INPUT ARGUEMENT VARIABLES
     a(.) = initial parameter set
     bl(.) = lower bound on parameters
     bu(.) = upper bound on parameters
     nopt = number of parameters to be optimized

  LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:
     ngs = number of complexes in the initial population
     npg = number of points in each complex
     npt = total number of points in initial population (npt=ngs*npg)
     nps = number of points in a sub-complex
     nspl = number of evolution steps allowed for each complex before
         complex shuffling
     mings = minimum number of complexes required, if the number of
         complexes is allowed to reduce as the optimization proceeds
     iseed = initial random seed
     iniflg = flag on whether to include the initial point in population
         = 0, not included
         = 1, included
     iprint = flag for controlling print-out after each shuffling loop
         = 0, print information on the best point of the population
         = 1, print information on every point of the population

  CONVERGENCE CHECK PARAMETERS
     maxn = max no. of trials allowed before optimization is terminated
     kstop = number of shuffling loops in which the criterion value must
         chang by the given percentage before optimization is terminated
     pcento = percentage by which the criterion value must change in
         given number of shuffling loops
     ipcnvg = flag indicating whether parameter convergence is reached
         (i.e., check if gnrng is less than 0.001)
         = 0, parameter convergence not satisfied
         = 1, parameter convergence satisfied

  LIST OF LOCAL VARIABLES
     x(.,.) = coordinates of points in the population
     xf(.) = function values of x(.,.)
     xx(.) = coordinates of a single point in x
     cx(.,.) = coordinates of points in a complex
     cf(.) = function values of cx(.,.)
     s(.,.) = coordinates of points in the current simplex
     sf(.) = function values of s(.,.)
     bestx(.) = best point at current shuffling loop
     bestf = function value of bestx(.)
     worstx(.) = worst point at current shuffling loop
     worstf = function value of worstx(.)
     xnstd(.) = standard deviation of parameters in the population
     gnrng = normalized geometric mean of parameter ranges
     lcs(.) = indices locating position of s(.,.) in x(.,.)
     bound(.) = bound on ith variable being optimized
     ngs1 = number of complexes in current population
     ngs2 = number of complexes in last population
     iseed1 = current random seed
     criter(.) = vector containing the best criterion values of the last
         10 shuffling loops
******************************************************************************/
void SCEUA::sceua
(
   double * a,
   double * bl,
   double * bu,
   int nopt,
   int maxn,
   int kstop,
   double pcento,
   int iseed,
   int ngs,
   int npg,
   int nps,
   int nspl,
   int mings,
   int iniflg,
   int iprint
)
{
   FILE * pOut = fopen("sce.out", "a");

   //LOCAL ARRAYS
   double ** x, * xx, * bestx, * worstx, * xf;
   double ** s, * sf, ** cx, * cf;
   double * xnstd, * bound, criter[20], * unit;
   int * lcs;   

   //allocate memory
   lcs = new int[nps];
   sf = new double[nps];
   xf = new double[ngs*npg];
   cf = new double[npg];
   x = new double * [ngs*npg];
   cx = new double * [npg];
   int i;
   for(i = 0; i < ngs*npg; i++)
   {
      x[i] = new double[nopt];
   }
   for(i = 0; i < npg; i++)
   {
      cx[i] = new double[nopt];
   }
   xx = new double[nopt];
   bestx = new double[nopt];
   worstx = new double[nopt];
   xnstd = new double[nopt];
   bound = new double[nopt];
   unit = new double[nopt];
   s = new double *[nps];
   for(i = 0; i < nps; i++)
   {
      s[i] = new double[nopt];
   }


   char **xname; 
   xname = new char*[nopt];
   for(i = 0; i < nopt; i++)
   {
     xname[i] = new char[50];   
     strncpy(xname[i], m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetName(), 9);
   }

   if(m_OutputMode != 2) printf("ENTER THE SCEUA SUBROUTINE --- \n");

   // INITIALIZE VARIABLES
   int nloop, loop, igs, nopt1, nopt2, itmp;
   nloop = 0;
   loop = 0;
   igs = 0;
   nopt1 = 8;
   if (nopt < 8) nopt1 = nopt;
   nopt2 = 12;
   if (nopt < 12) nopt2 = nopt;

   //INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER
   int iseed1;
   iseed1 = -abs(iseed);

   //COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
   int npt, ngs1, npt1;
   npt = ngs * npg; 
   ngs1 = ngs; 
   npt1 = npt; 

   fprintf(pOut, "\
  ==================================================\n\
  ENTER THE SHUFFLED COMPLEX EVOLUTION GLOBAL SEARCH\n\
  ==================================================\n\n\n");

   if(m_OutputMode != 2) printf(" ***  Evolution Loop Number %d\n", nloop);

   //COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
   int j;
   for(j = 1; j <= nopt; j++)
   {
      bound[j-1] = bu[j-1] - bl[j-1];
      unit[j-1] = 1.0;
   }

   //COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT   
   double fa;
   //handle warm start
   if(m_pModel->CheckWarmStart() == true)
   {
      WarmStart();
      m_pModel->GetParamGroupPtr()->ReadParams(a);
   }
   // handle parameter extraction
   if(m_pModel->GetParamGroupPtr()->CheckExtraction() == true)
   {
      m_pModel->GetParamGroupPtr()->ReadParams(a);
   }
   m_pModel->GetParamGroupPtr()->WriteParams(a);
   fa = m_pModel->Execute();
   if (fa < m_fSaved)
   {
      m_fSaved = fa;
      m_pModel->SaveBest(0);
   }

   //write initial config.
   int nleft = m_Budget - m_pModel->GetCounter();

   double eb = (double)(m_pModel->GetCounter())/(double)m_Budget; //elapsed budget
   
   m_pStatus.curIter = 0;
   m_pStatus.maxIter = m_Budget;
   m_pStatus.pct = (float)100.00*((float)1.00-(float)nleft/(float)m_Budget);
   m_pStatus.numRuns = m_pModel->GetCounter();
   m_pModel->GetParamGroupPtr()->WriteParams(m_pParams);
   WriteRecord(m_pModel, 0, fa, m_pStatus.pct);
   m_CurIter++;
   WriteStatus(&m_pStatus);

   //PRINT THE INITIAL POINT AND ITS CRITERION VALUE
   fprintf(pOut, "\
*** PRINT THE INITIAL POINT AND ITS CRITERION VALUE ***\n\n\
 CRITERION    ");
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut, "%-9s    ", xname[i]);
   }
   fprintf(pOut, "\n  %8.0lf     ", fa);
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut, "%5.3lf     ", a[i]);
   }
   fprintf(pOut, "\n\n\n");
   
   // GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
   if (iniflg == 1)
   {
      for(j = 0; j < nopt; j++)
         x[0][j] = a[j];
      xf[0] = fa;
      WriteInnerEval(WRITE_SCE, npt, '.');
      WriteInnerEval(1, npt, '.');
   }/* end if() */
   // ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO x(1,.)
   else
   {
      getpnt(nopt,1,&iseed1,xx,bl,bu,unit,bl);
      eb = (double)(m_pModel->GetCounter())/(double)m_Budget;
      for(j = 0; j < nopt; j++)
      {
         xx[j] = TelescopicCorrection(bl[j], bu[j], bestx[j], eb, xx[j]);
         x[0][j] = xx[j];
      }
      WriteInnerEval(WRITE_SCE, npt, '.');
      WriteInnerEval(1, npt, '.');
      m_pModel->GetParamGroupPtr()->WriteParams(xx);
      m_pModel->PerformParameterCorrections();
      for(j = 0; j < nopt; j++)
      {
         xx[j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal();
         x[0][j] = xx[j];
      }
      xf[0] = m_pModel->Execute();
      if (xf[0] < m_fSaved)
      {
         m_fSaved = xf[0];
         m_pModel->SaveBest(0);
      }
   }

   //use initial point if it's better than the randomly generated one
   if(fa < xf[0])
   {
      fprintf(pOut, "THE INITIAL POINT IS BETTER THAN THE RANDOM STARTING POINT. USING IT INSTEAD.");

      for(j = 0; j < nopt; j++)
         x[0][j] = a[j];
      xf[0] = fa;
   }

   int icall;
   icall = 1;
   if (icall >= maxn) goto label_9000;

   //GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
   //SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES
label_restart:
   for(i = 1; i < npt1; i++)
   {
      getpnt(nopt,1,&iseed1,xx,bl,bu,unit,bl);
      eb = (double)(m_pModel->GetCounter())/(double)m_Budget;
      for(j = 0; j < nopt; j++)
      {
         xx[j] = TelescopicCorrection(bl[j], bu[j], bestx[j], eb, xx[j]);
         x[i][j] = xx[j];
      }
      WriteInnerEval(i+1, npt, '.');
      m_pModel->GetParamGroupPtr()->WriteParams(xx);
      m_pModel->PerformParameterCorrections();
      for(j = 0; j < nopt; j++)
      {
         xx[j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal();
         x[i][j] = xx[j];
      }
      xf[i] = m_pModel->Execute();
      if (xf[i] < m_fSaved)
      {
         m_fSaved = xf[i];
         m_pModel->SaveBest(0); 
      }

      icall = icall + 1;
      if (icall >= maxn)
      {
         break;
      }
   }/* end for() */
   WriteInnerEval(WRITE_ENDED, npt, '.');

   // ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
   sort(npt1,nopt,x,xf);

   for(j = 0; j < nopt; j++)
   {
      bestx[j] = x[0][j];
      worstx[j] = x[npt1-1][j];
   }
   double bestf, worstf;
   bestf = xf[0];
   worstf = xf[npt1-1];

   // COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
   double gnrng;
   int ipcnvg;
   parstt(npt1,nopt,x,xnstd,bound,&gnrng,&ipcnvg);

   // PRINT THE RESULTS FOR THE INITIAL POPULATION
   fprintf(pOut,"\
**** PRINT THE RESULTS OF THE SCE SEARCH ***\n\n\
 LOOP TRIALS COMPLXS  BEST F   WORST F   PAR RNG         ");
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut,"%-9s ", xname[i]);
   }
   fprintf(pOut,"\n");
   fprintf(pOut," %4d %6d %7d  %6.2lf  %9.3E  %8.3lf      ",
           nloop,icall,ngs1,bestf,worstf,gnrng);
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut,"%6.3lf    ", bestx[i]);
   }
   fprintf(pOut,"\n");

   if (iprint == 1)
   {
      fprintf(pOut, "POPULATION AT LOOP (%d)\n", nloop);
      for(i = 0; i < npt1; i++)
      {
         fprintf(pOut, "%8.3lf    ", xf[i]);
         for(j = 0; j < nopt; j++)
         {
            fprintf(pOut, "%8.3lf    ", x[i][j]);
         }
         fprintf(pOut, "\n");
      }/* end for() */
   }/* end if() */

   if (icall >= maxn) goto label_9000;
   if (ipcnvg == 1) goto label_9200;

   // BEGIN THE MAIN LOOP ----------------
 label_1000:

   nleft = m_Budget - m_pModel->GetCounter();
   m_pStatus.curIter = nloop+1;
   if(IsQuit() == true){ goto label_9999;}
   if(nleft <= 0){ m_pStatus.pct = 100.00;  goto label_9000;}

   nloop = nloop + 1;

   if(m_OutputMode != 2) printf(" ***  Evolution Loop Number %d\n",nloop); 
   
   //BEGIN LOOP ON COMPLEXES
   for(igs = 1; igs <= ngs1; igs++)
   {
      // ASSIGN POINTS INTO COMPLEXES
      int k1, k2;
      for(k1 = 1; k1 <= npg; k1++) 
      {
        k2 = (k1-1) * ngs1 + igs; 
        for(j = 1; j <= nopt; j++)
        {
            cx[k1-1][j-1] = x[k2-1][j-1]; 
        } //end for()
        cf[k1-1] = xf[k2-1];
      } //end for()
      // BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
      int tmp = 0;
      WriteInnerEval(WRITE_SCE, m_NumEvoSteps, '.');

      for(loop = 0; loop < nspl; loop++) 
      {
         // CHOOSE A SUB-COMPLEX (nps points) ACCORDING TO A LINEAR
         // PROBABILITY DISTRIBUTION
         if (nps == npg)
         {
            int k;
            for(k = 0; k < nps; k++)
            {
               lcs[k] = k;
            } // end do
            goto label_85; 
         } // end if

         double myrand;
         myrand = UniformRandom();
         itmp = (int)(npg + 0.5 - sqrt(pow((npg+0.5),2.00) - npg*(npg+1.00)*myrand));
         if(itmp >= npg) itmp = npg-1;
         lcs[0] = itmp;

         int k, lpos;
         for(k = 2; k <= nps; k++) 
         {
label_60:
            myrand = UniformRandom(); 
            lpos = (int)(npg + 0.5 - sqrt(pow((npg+0.5),2.00) - npg*(npg+1.00)*myrand));
            if(lpos >= npg) lpos = npg-1;

            for(k1 = 1; k1 <= k-1; k1++) 
            {
               if (lpos == lcs[k1-1]) goto label_60; 
            } // end do
            lcs[k-1] = lpos; 
         } // end do

         // ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
         sort(nps,lcs);

         // CREATE THE SUB-COMPLEX ARRAYS
label_85: 
         for(k = 1; k <= nps; k++)
         {
            for(j = 1; j <= nopt; j++) 
            {
               s[k-1][j-1] = cx[lcs[k-1]][j-1]; 
            } // end do
            sf[k-1] = cf[lcs[k-1]];
         } // end do

         // USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
         cce(nopt,nps,s,sf,bl,bu,xnstd,&tmp,maxn,&iseed1);

         // IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
         // INTO THE COMPLEX
         for(k = 1; k <= nps; k++) 
         {
            for(j = 1; j <= nopt; j++) 
            {
               cx[lcs[k-1]][j-1] = s[k-1][j-1]; 
            } // end do
            cf[lcs[k-1]] = sf[k-1]; 
         } //end do

         // SORT THE POINTS
         sort(npg,nopt,cx,cf);

         //IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
         if (icall >= maxn) break; 
         // END OF INNER LOOP ------------
      } /* end for() */

      WriteInnerEval(WRITE_ENDED, m_NumEvoSteps, '.');
      icall += tmp;

      // REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
      for(k1 = 1; k1 <= npg; k1++)
      {
         k2 = (k1-1) * ngs1 + igs; 
         for(j = 1; j <= nopt; j++) 
         {
            x[k2-1][j-1] = cx[k1-1][j-1]; 
         } // end do
         xf[k2-1] = cf[k1-1]; 
      } // end do
      if (icall >= maxn) break; 
      //END LOOP ON COMPLEXES
   } /* end for() */ 

   // RE-SORT THE POINTS
   sort(npt1,nopt,x,xf);

   // RECORD THE BEST AND WORST POINTS
   for(j = 0; j < nopt; j++) 
   {
      m_pParams[j] = bestx[j] = x[0][j];
      worstx[j] = x[npt1-1][j]; 
   } //end do
   m_Best = bestf = xf[0];
   worstf = xf[npt1-1];

   // TEST THE POPULATION FOR PARAMETER CONVERGENCE
   parstt(npt1,nopt,x,xnstd,bound,&gnrng,&ipcnvg);

   // PRINT THE RESULTS FOR CURRENT POPULATION
   m_pModel->GetParamGroupPtr()->WriteParams(m_pParams);
   nleft = m_Budget - m_pModel->GetCounter();
   m_pStatus.pct = (float)100.00*((float)1.00-(float)nleft/(float)m_Budget);
   m_pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&m_pStatus);
   WriteRecord(m_pModel, nloop, m_Best, m_pStatus.pct);
   m_CurIter++;

   if((nloop%5) == 0)
   {
      fprintf(pOut,"\
 LOOP TRIALS COMPLXS  BEST F   WORST F   PAR RNG         ");
      for(i = 0; i < nopt; i++)
      {
         fprintf(pOut,"%-9s ", xname[i]);
      }
      fprintf(pOut,"\n");
   }
   fprintf(pOut," %4d %6d %7d  %6.2lf  %9.3E  %8.3lf      ",
           nloop,icall,ngs1,bestf,worstf,gnrng);
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut,"%6.3lf    ", bestx[i]);
   }
   fprintf(pOut,"\n");

   if (iprint == 1)
   {
      fprintf(pOut, "POPULATION AT LOOP (%d)\n", nloop);
      for(i = 0; i < npt1; i++)
      {
         fprintf(pOut, "%8.3lf    ", xf[i]);
         for(j = 0; j < nopt; j++)
         {
            fprintf(pOut, "%8.3lf    ", x[i][j]);
         }
         fprintf(pOut, "\n");
      }/* end for() */
   }/* end if() */

   // TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
   if (icall >= maxn) goto label_9000; 

   // COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
   criter[19] = bestf; 
   if (nloop < (kstop+1)) goto label_132; 
   double denomi, timeou;
   denomi = fabs(criter[19-kstop] + criter[19]) / 2.0; 
   timeou = fabs(criter[19-kstop] - criter[19]) / denomi; 
   if (timeou < pcento) goto label_9100; 
label_132: 
   int l;
   for(l=0; l < 19; l++) 
   {
        criter[l] = criter[l+1];
   } //end do


   //IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
   if (ipcnvg == 1) goto label_9200; 

   //NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH

   //CHECK FOR COMPLEX NUMBER REDUCTION
   int ngs2;
   if (ngs1 > mings) 
   {
        ngs2 = ngs1; 
        ngs1 = ngs1 - 1; 
        npt1 = ngs1 * npg; 
        comp(nopt,npt1,ngs1,ngs2,npg,x,xf,cx,cf); 
   } //end if

   // END OF MAIN LOOP -----------
   goto label_1000;

   // SEARCH TERMINATED
label_9000: 

      fprintf(pOut, "\
*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE\n\
    LIMIT ON THE MAXIMUM NUMBER OF TRIALS (%d)\n\
    WAS EXCEEDED.  SEARCH WAS STOPPED AT %d SUB-COMPLEX\n\
    OF COMPLEX %d IN SHUFFLING LOOP %d ***\n\n",
    maxn, loop, igs, nloop);

      goto label_9999; 

label_9100:

      fprintf(pOut, "\
*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION\n\
    VALUE HAS NOT CHANGED %5.2lf PERCENT IN %d\n\
    SHUFFLING LOOPS ***\n\n", pcento*100.0,kstop);

      goto label_9999; 

 label_9200:
      fprintf(pOut, "\
 *** OPTIMIZATION RESTARTED BECAUSE THE POPULATION HAS\n\
     CONVERGED INTO %5.2lf PERCENT OF THE FEASIBLE SPACE ***\n\n",
      gnrng*100.0); 
      goto label_restart;

 label_9999:
      // PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      fprintf(pOut, "\
*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS CRITERION VALUE ***\n\n\
 CRITERION        ");

      for(i = 0; i < nopt; i++)
      {
         fprintf(pOut,"%-9s ", xname[i]);
      }
      fprintf(pOut,"\n%6.3lf    ", bestf);
      for(i = 0; i < nopt; i++)
      {
         fprintf(pOut,"%6.3lf    ", bestx[i]);
      }
      fprintf(pOut,"\n");

   fclose(pOut);
   /* clean up memory */
   for(i = 0; i < nopt; i++)
   {
     delete [] xname[i];   
   }
   delete [] xname;
   for(i = 0; i < ngs*npg; i++)
   {
      delete [] x[i];
   }
   for(i = 0; i < npg; i++)
   {
      delete [] cx[i];
   }
   delete [] x;
   delete [] cx;
   delete [] xx;
   delete [] bestx;
   delete [] worstx;
   delete [] xnstd;
   delete [] bound;
   delete [] unit;
   for(i = 0; i < nps; i++)
   {
      delete [] s[i];
   }
   delete [] s;
   delete [] sf;
   delete [] lcs;
   delete [] cf;
   delete [] xf;
} /* end sceua() */

/******************************************************************************
cce()

ALGORITHM TO GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
******************************************************************************/
void SCEUA::cce
(
   int nopt,
   int nps,
   double ** s,
   double * sf,
   double * bl,
   double * bu,
   double * xnstd,
   int * icall,
   double maxn,
   int * iseed
)
{
   //SUB-COMPLEX VARIABLES
   const double c1=0.8;
   const double c2=0.4;

   /* ----------------------------------------------------
   LIST OF LOCAL VARIABLES
      sb(.) = the best point of the simplex
      sw(.) = the worst point of the simplex
      w2(.) = the second worst point of the simplex
      fw = function value of the worst point
      ce(.) = the centroid of the simplex excluding wo
      snew(.) = new point generated from the simplex
      iviol = flag indicating if constraints are violated
            = 1 , yes
            = 0 , no
   ----------------------------------------------------- */
   double * sw, * sb, *ce, *snew;
   sw = new double[nopt];
   sb = new double[nopt];
   ce = new double[nopt];
   snew = new double[nopt];

   //EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
   int n = nps;
   int m = nopt;
   double alpha = 1.0;
   double beta = 0.5;

   /* ---------------------------------------------------
   IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
   COMPUTE THE CENTROID ce OF THE REMAINING POINTS
   COMPUTE step, THE VECTOR BETWEEN wo AND ce
   IDENTIFY THE WORST FUNCTION VALUE fw
   --------------------------------------------------- */
   int i, j;
   for(j = 0; j < m; j++)
   {
      sb[j] = s[0][j]; 
      sw[j] = s[n-1][j]; 
      ce[j] = 0.0; 
      for(i = 0; i < n-1; i++) 
      {
         ce[j] = ce[j] + s[i][j]; 
      } //end do
     ce[j] = ce[j]/(double)(n-1);
   } //end do
   double fw;
   fw = sf[n-1]; 

   //COMPUTE THE NEW POINT snew
   //FIRST TRY A REFLECTION STEP
   for(j = 0; j < m; j++) 
   {
      snew[j] = ce[j] + alpha * (ce[j] - sw[j]);
   } //end do

   //CHECK IF snew SATISFIES ALL CONSTRAINTS
   int ibound;
   chkcst(nopt,snew,bl,bu,&ibound); 

   /* ------------------------------------------------------------------
   snew IS OUTSIDE THE BOUND,
   CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
   A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
   AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
   ------------------------------------------------------------------ */
   if (ibound >= 1) getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb);

   //COMPUTE THE FUNCTION VALUE AT snew
   double fnew;
   WriteInnerEval(*icall+1, m_NumEvoSteps, '.');
   double eb = (double)(m_pModel->GetCounter())/(double)m_Budget;
   for(j = 0; j < m; j++) snew[j] = TelescopicCorrection(bl[j], bu[j], sb[j], eb, snew[j]);
   m_pModel->GetParamGroupPtr()->WriteParams(snew);
   m_pModel->PerformParameterCorrections();
   for(j = 0; j < m; j++) snew[j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal();
   fnew = m_pModel->Execute(); 
   if (fnew < m_fSaved)
   {
      m_fSaved = fnew;
      m_pModel->SaveBest(0);
   }

   *icall = *icall + 1;

   //COMPARE fnew WITH THE WORST FUNCTION VALUE fw
   //fnew IS LESS THAN fw, ACCEPT THE NEW POINT snew AND RETURN
   if (fnew <= fw) goto label_2000;
   if (*icall >= maxn) goto label_3000;

   //fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
   for(j = 0; j < m; j++)
   {
      snew[j] = ce[j] - beta * (ce[j] - sw[j]);
   } //end do

   //COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
   eb = (double)(m_pModel->GetCounter())/(double)m_Budget;
   for(j = 0; j < m; j++) snew[j] = TelescopicCorrection(bl[j], bu[j], sb[j], eb, snew[j]);
   WriteInnerEval(*icall+1, m_NumEvoSteps, '.');
   m_pModel->GetParamGroupPtr()->WriteParams(snew);
   m_pModel->PerformParameterCorrections();
   for(j = 0; j < m; j++) snew[j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal();
   fnew = m_pModel->Execute();
   if (fnew < m_fSaved)
   {
      m_fSaved = fnew;
      m_pModel->SaveBest(0);
   }

   *icall = *icall + 1;

   //COMPARE fnew TO THE WORST VALUE fw
   //IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN
   if (fnew <= fw) goto label_2000; 
   if (*icall >= maxn) goto label_3000; 

   /* ---------------------------------------------------------------------
   IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
   ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
   AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
   --------------------------------------------------------------------- */
   getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb);

   //COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
   eb = (double)(m_pModel->GetCounter())/(double)m_Budget;
   for(j = 0; j < m; j++) snew[j] = TelescopicCorrection(bl[j], bu[j], sb[j], eb, snew[j]);
   WriteInnerEval(*icall+1, m_NumEvoSteps, '.');
   m_pModel->GetParamGroupPtr()->WriteParams(snew);
   m_pModel->PerformParameterCorrections();
   for(j = 0; j < m; j++) snew[j] = m_pModel->GetParamGroupPtr()->GetParamPtr(j)->GetEstVal();
   fnew = m_pModel->Execute();
   if (fnew < m_fSaved)
   {
      m_fSaved = fnew;
      m_pModel->SaveBest(0);
   }

   *icall = *icall + 1;

   //REPLACE THE WORST POINT BY THE NEW POINT
label_2000:
   for(j = 0; j < m; j++)
   {
      s[n-1][j] = snew[j];
   } //end do
   sf[n-1] = fnew;

label_3000:
   //free up memory  
   delete [] sw;
   delete [] sb;
   delete [] ce;
   delete [] snew;
} /* end cce() */

/******************************************************************************
getpnt()

This subroutine generates a new point within feasible region

x(.) = new point
xi(.) = focal point
bl(.) = lower bound
bu(.) = upper bound
std(.) = standard deviation of probability distribution
idist = probability flag
      = 1 - uniform distribution
      = 2 - Gaussian distribution
******************************************************************************/
void SCEUA::getpnt
(
   int nopt,
   int idist,
   int * iseed,
   double * x,
   double * bl,
   double * bu,
   double * std,
   double * xi
)
{
   int ibound;
   int j;
   double myrand;

   int icount = m_pModel->GetCounter();
label_1:
   for(j = 0; j < nopt; j++) //1   do j=1, nopt
   {
label_2:
      if (idist == 1)
      {
         myrand = UniformRandom();
      }
      else if (idist == 2)
      {
         myrand = GaussRandom();
      }
      else
      {
         printf("unknown distribution!\n");
      }

      x[j] = xi[j] + std[j] * myrand * (bu[j] - bl[j]);

      FILE * pFile = fopen("getpnt.txt", "a+");
      fprintf(pFile, "%d\tx[%d]:%E\txi[%d]:%e\tstd[%d]:%E\tmyrand : %E\tbu[%d]:%E\tbl[%d]:%E\n",
                     icount, j+1, x[j], j+1, xi[j], j+1, std[j], myrand, j+1, bu[j], j+1, bl[j]);
      fclose(pFile);

      //Check explicit constraints
      chkcst(1,&x[j],&bl[j],&bu[j],&ibound);
      if (ibound >= 1)
      {
         goto label_2;
      }
   } //end do

   //Check implicit constraints    
   chkcst(nopt,x,bl,bu,&ibound); 
   if (ibound >= 1)
   {
      goto label_1; 
   }
}/* end getpnt() */

/******************************************************************************
parstt()

SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
******************************************************************************/
void SCEUA::parstt
(
   int npt,
   int nopt,
   double ** x,
   double * xnstd,
   double * bound,
   double * gnrng,
   int * ipcnvg
)
{
   double * xmax, * xmin, * xmean;
   xmax = new double[nopt];
   xmin = new double[nopt];
   xmean = new double[nopt];

   const double delta = 1.0e-20;
   const double peps = m_Peps;

   //COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
   double gsum, xsum1, xsum2;
   int i,k;
   gsum = 0.0; 
   for(k = 0; k < nopt; k++) 
   {
      xmax[k] = -1.0E+20; 
      xmin[k] = 1.0E+20;
      xsum1 = 0.0; 
      xsum2 = 0.0; 

      for(i = 0; i < npt; i++) 
      {
         xmax[k] = MyMax(x[i][k], xmax[k]); 
         xmin[k] = MyMin(x[i][k], xmin[k]); 
         xsum1 = xsum1 + x[i][k]; 
         xsum2 = xsum2 + x[i][k]*x[i][k];
      } //end do

      xmean[k] = xsum1/(double)npt; 
      xnstd[k] = (xsum2 / (double)npt - xmean[k]*xmean[k]);

      if (xnstd[k] <= delta) xnstd[k] = delta; 
      xnstd[k] = sqrt(xnstd[k]);
      xnstd[k] = xnstd[k] / bound[k];

      gsum = gsum + log(delta + (xmax[k]-xmin[k])/bound[k]);
   } //end do
   *gnrng = exp(gsum/(double)nopt);

   //CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
   *ipcnvg = 0;
   if (*gnrng <= peps) 
   {
      *ipcnvg = 1;
   } //end if

   delete [] xmax;
   delete [] xmin;
   delete [] xmean;

}/* end parstt() */

/******************************************************************************
comp()

THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
******************************************************************************/
void SCEUA::comp
(
   int n,
   int npt,
   int ngs1,
   int ngs2,
   int npg,
   double ** a,
   double * af,
   double ** b,
   double * bf
)
{
   int i, igs, ipg, k1,  k2;
   for(igs = 1; igs <= ngs1; igs++)
   {
      for(ipg = 1; ipg <= npg; ipg++)
      {
         k1=(ipg-1)*ngs2 + igs; 
         k2=(ipg-1)*ngs1 + igs;
         while(k2 > npg) k2 -= npg;
         for(i = 1; i <= n; i++)
         {
            b[k2-1][i-1] = a[k1-1][i-1];
         } //end do
         bf[k2-1] = af[k1-1];
      } //end do
   } //end do

   int j, jj;
   for(j = 0; j < npt; j++) 
   {
      jj = j % npg;
      for(i = 0; i < n; i++)
      {
         a[j][i] = b[jj][i]; 
      } //end do
      af[j] = bf[jj];
   } //end do

} /* end comp() */

/******************************************************************************
sort()

Sort a complex of parameter sets in ascending order of cost function.
******************************************************************************/
void SCEUA::sort
(
   int n,
   int m,
   double ** rb,
   double * ra
)
{
   //n = number of paramter sets
   //m = number of parameters
   //rb = list of parameter sets
   //ra = list of costs

   for(int p = 0; p < n; p++)
   {
      double fp = ra[p];
      for(int i = p + 1; i < n; i++)
      {
         double fi = ra[i];
         //swap if ith entry is better than pth
         if(fi < fp)
         {
            ra[p] = fi;
            ra[i] = fp;
            for(int j = 0; j < m; j++)
            {
               double t = rb[i][j];
               rb[i][j] = rb[p][j];
               rb[p][j] = t;
            }/* end for() */
         }/* end if() */
      }/* end for() */
   }/* end for() */
} /* end sort() */

/******************************************************************************
sort()

Sort a list of integers in ascending order.
******************************************************************************/
void SCEUA::sort
(
   int n,
   int * ra
)
{
   for(int p = 0; p < n; p++)
   {
      int fp = ra[p];
      for(int i = p + 1; i < n; i++)
      {
         int fi = ra[i];
         //swap if ith entry is better than pth
         if(fi < fp)
         {
            ra[p] = fi;
            ra[i] = fp;
         }/* end if() */
      }/* end for() */
   }/* end for() */
} /* end sort() */

/******************************************************************************
chkcst()

     This subroutine check if the trial point satisfies all
     constraints.

     ibound - violation indicator
            = -1 initial value
            = 0  no violation
            = 1  violation
     nopt = number of optimizing variables
     ii = the ii'th variable of the arrays x, bl, and bu
******************************************************************************/
void SCEUA::chkcst
(
   int nopt,
   double * x,
   double * bl,
   double * bu,
   int * ibound
)
{
   *ibound = -1;

   //Check if explicit constraints are violated
   for(int ii=1; ii<=nopt; ii++)
   {
      if ((x[ii-1] < bl[ii-1]) || (x[ii-1] > bu[ii-1])) goto label10;
   }
   if (nopt == 1) goto label9;


//     Check if implicit constraints are violated
//     (no implicit constraints for this function)
//
//     No constraints are violated
//      
label9:    *ibound = 0;
      return;

//    At least one of the constraints are violated
label10:   *ibound = 1;
      return;
}/* end chkcst() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename, then write the 
configuration info. to the file "sce.in" (maintains compatibility with SCE
fortran implementation).
******************************************************************************/
void SCEUA::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   int i;
   char * line;
   char tmp[DEF_STR_SZ], tmp2[DEF_STR_SZ];

   //assign defaults
   m_np = m_pModel->GetParamGroupPtr()->GetNumParams();
   m_Budget = 10000; //MAXN
   m_Kstop = 5; //KSTOP
   m_Pcento = 0.01; //PCENTO
   m_Peps = 1.0E-3; //peps
   m_NumComplexes = (int)(sqrt((double)m_np)); //NGS
   m_Seed = 1969; //ISEED
   m_UserConfig = 1; //IDEFLT
   m_PtsPerComplex = 2*m_np + 1; //NPG
   m_PtsPerSubComplex = m_np+1; //NPS
   m_NumEvoSteps = m_PtsPerComplex; //NSPL
   m_MinComplexes = m_NumComplexes; //MINGS
   m_UseInitPt = 1; //INIFLG
   m_OutputMode = 2; //IPRINT
   m_bUseInitPt = false; //INIFLG

   //allocate initial parameter configuration
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();
   m_np = pGroup->GetNumParams();
   NEW_PRINT("double", m_np);
   m_pParams = new double[m_np];
   MEM_CHECK(m_pParams);

   NEW_PRINT("double", m_np);
   m_pLower = new double[m_np];
   MEM_CHECK(m_pParams);

   NEW_PRINT("double", m_np);
   m_pUpper = new double[m_np];
   MEM_CHECK(m_pParams);

   for(i = 0; i < m_np; i++)
   {
      m_pParams[i] = pGroup->GetParamPtr(i)->GetEstVal();
      m_pLower[i] = pGroup->GetParamPtr(i)->GetLwrBnd();
      m_pUpper[i] = pGroup->GetParamPtr(i)->GetUprBnd();
   }/* end for() */

   //read in SCEUA configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open SCEUA config. file. Using Defaults");      
      return;
   }/* end if() */   

   if(CheckToken(pFile, "RandomSeed", GetInFileName()) == true)
   {
      fclose(pFile);
      m_Seed = GetRandomSeed();
      pFile = fopen(pFileName, "r");
   }
   rewind(pFile);

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginSCEUA", pFileName) == true)
   {
      FindToken(pFile, "EndSCEUA", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginSCEUA", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndSCEUA") == NULL)
      {         
         if(strstr(line, "Budget") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_Budget); 
            if(m_Budget < 100)
            {
               LogError(ERR_FILE_IO, "Invalid SCEUA budget. Defaulting to 100.");
               m_Budget = 100;
            }
         }/*end if() */
         else if(strstr(line, "LoopStagnationCriteria") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_Kstop); 
         }
         else if(strstr(line, "PctChangeCriteria") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Pcento); 
         }
         else if(strstr(line, "PopConvCriteria") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Peps); 
         }          
         else if(strstr(line, "NumComplexes") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumComplexes); 
         }
         else if(strstr(line, "NumPointsPerComplex") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_PtsPerComplex); 
         }
         else if(strstr(line, "NumPointsPerSubComplex") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_PtsPerSubComplex); 
         }
         else if(strstr(line, "NumEvolutionSteps") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumEvoSteps); 
         }
         else if(strstr(line, "MinNumOfComplexes") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MinComplexes); 
         }
         else if(strstr(line, "UseInitialPoint") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2); 
            MyStrLwr(tmp2);
            if(strcmp(tmp2, "yes") == 0) m_bUseInitPt = true;
         }

         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   
   fclose(pFile);

   //create sce.in file
   int iniflg = 0;
   if(m_bUseInitPt) iniflg = 1;
   FILE * pOut = fopen("sce.in", "w");
   fprintf(pOut, "%d  %d  %lf  %d  %d  1\n", 
           m_Budget, m_Kstop, m_Pcento, m_NumComplexes, m_Seed);
   fprintf(pOut, "%d  %d  %d  %d  %d  2\n", 
           m_PtsPerComplex, m_PtsPerSubComplex, m_NumEvoSteps, m_MinComplexes, iniflg);
   for(i = 0; i < m_np; i++)
   {
      fprintf(pOut, "%E %E %E\n", m_pParams[i], m_pLower[i], m_pUpper[i]);
   }
   fclose(pOut);
} /* end InitFromFile() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void SCEUA::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm                : Shuffled Complex Evolution (SCE)\n");
   fprintf(pFile, "Budget                   : %d\n", m_Budget);
   fprintf(pFile, "Loop Stagnation Criteria : %d\n", m_Kstop); 
   fprintf(pFile, "Pct Change Criteria      : %lf\n", m_Pcento); 
   fprintf(pFile, "Number of Complexes      : %d\n", m_NumComplexes); 
   fprintf(pFile, "Points Per Complex       : %d\n", m_PtsPerComplex); 
   fprintf(pFile, "Points Per Sub-Complex   : %d\n", m_PtsPerSubComplex); 
   fprintf(pFile, "Num. of Evolution Steps  : %d\n", m_NumEvoSteps); 
   fprintf(pFile, "Min. Num. of Complexes   : %d\n", m_MinComplexes); 
  
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
SCEUA_Program()

Calibrate or optimize the model using SCE.
******************************************************************************/
void SCEUA_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("SCEUA", 1);
   SCEUA * SCE = new SCEUA(model);
   MEM_CHECK(SCE);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { SCE->Calibrate(); }
   else { SCE->Optimize(); }

   delete SCE;
   delete model;
} /* end SCEUA_Program() */
