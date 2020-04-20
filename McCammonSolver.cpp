/******************************************************************************
File     : McCammonSolver.cpp
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The McCammonSolver class solves the non-linear McCammon equation for handling
errors in both variables (q and C).

See McCammon.cpp for details.

Version History
07-28-07    lsm   added copyright information and initial comments.
******************************************************************************/
#include <string.h>

#include "McCammonSolver.h"
#include "ObservationGroup.h"
#include "Observation.h"

#include "Isotherms.h"
#include "IsoParse.h"
#include "Utility.h"
#include "Exception.h"

double * gC0 = NULL;

/******************************************************************************
CTOR

Constructs a McCammonSolver class using the given Isotherm.
******************************************************************************/
McCammonSolver::McCammonSolver(IsothermABC * pIso)
{
   int i;
   m_pWc = NULL;
   m_pq  = NULL;
   m_pWq = NULL;
   m_MaxIters = 50;
   
   m_pIso = pIso;

   //share certain data members with the Isotherm class
   m_pC = pIso->GetPtrToC(&m_NumOut);
   m_pOutFile = pIso->GetPtrToOutFile();

   gC0 = new double[m_NumOut];

   // determine min and max concentrations (these will bound the search)
   m_Cupr = m_Clwr = m_pC[0];
   for(i = 0; i < m_NumOut; i++)
   {
      if(m_pC[i] < m_Clwr)
      {
         m_Clwr = m_pC[i];
      }
      if(m_pC[i] > m_Cupr)
      {
         m_Cupr = m_pC[i];
      }
      gC0[i] = m_pC[i]; //save initial concentrations (for diskless mode)
   }
   /* ----------------------------------------------
   Preferred lower bound is a concentration of zero.
   But this can cause divide-by-zero errors for 
   certain isotherms. So do some testing here to see
   if concentration of 0 is ok, otherwise use 1/10th
   of the lowest observed concentration value.   
   ---------------------------------------------- */
   double tmp1, tmp2;
   tmp1 = m_pIso->q(0);
   tmp2 = m_pIso->dqdc(0);
   if((CheckOverflow(tmp1) == true) || (CheckOverflow(tmp2) == true))
   {
      m_Clwr /= 10.00;
      if(m_Clwr < 1E-10) m_Clwr = 1E-10;
      //printf("Orear() Overflow condition detected and corrected");
   }
   else
   {
      m_Clwr = 0.00;
   }
   m_Cupr *= 2.00;

   IncCtorCount();
} /* end McCammonSolver::CTOR */

/******************************************************************************
DTOR

Destroys a McCammonSolver class.
******************************************************************************/
void McCammonSolver::Destroy(void)
{   
   delete [] m_pWc;
   delete [] m_pq;
   delete [] m_pWq;
   delete [] gC0;
   IncDtorCount();
} /* end McCammonSolver::DTOR */

/******************************************************************************
McCammonSolver::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void McCammonSolver::Compute(void)
{
   FILE * pFile;
   int i;

   /* -----------------------------------------------------------------
   For each data point, Calculate C that minimizes the nonlinear 
   equation, as defined in F(). Since there is  only one equation for 
   each data point, we use a simple "brute-force" bisection method.
   ----------------------------------------------------------------- */
   // perform bisection search
   for(i = 0; i < m_NumOut; i++)
   {
      m_pC[i] = BisectionSearch(i);
   }/* end for() */

   //Utilize isotherm to write out optimal C, q
   m_pIso->Compute();

   //write out Orear settings (Ct, S, and V)
   pFile = fopen(m_pOutFile, "a");
   fprintf(pFile, "\nSolutionMethod McCammon\n");
   fprintf(pFile, "Max Bisections %d\n", m_MaxIters);

   fprintf(pFile, "i   Aqueous Weight  Sorbed Weight\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%02d  %E   %E\n",i,m_pWc[i], m_pWq[i]);
   }   
   fclose(pFile);
}/* end Compute() */

/******************************************************************************
McCammonSolver::Compute()

Compute output vales and write them to the ObservationGroup.
******************************************************************************/
void McCammonSolver::Compute(ObservationGroup * pObs)
{
   int i;

   for(i = 0; i < m_NumOut; i++)
   {
      m_pC[i] = gC0[i]; //restore measured concentrations
   }

   /* -----------------------------------------------------------------
   For each data point, Calculate C that minimizes the nonlinear 
   equation, as defined in F(). Since there is  only one equation for 
   each data point, we use a simple "brute-force" bisection method.
   ----------------------------------------------------------------- */
   // perform bisection search
   for(i = 0; i < m_NumOut; i++)
   {
      m_pC[i] = BisectionSearch(i);
   }/* end for() */

   //Utilize isotherm to set simulated q
   m_pIso->Compute(pObs);

   //set simulated C
   for(i = 0; i < m_NumOut; i++){ pObs->GetObsPtr(i+m_NumOut)->SetComputedVal(m_pC[i]); }
}/* end Compute() */

/******************************************************************************
McCammonSolver::BisectionSearch()

Simple bisection search. Evaluates at least two new points each iteration and 
reduces search space by 50%.

Returns optimial concentration for the i-th observation point.
******************************************************************************/
double McCammonSolver::BisectionSearch(int i)
{
   double qobs, Cobs, wq, wc;
   double Cupr, Clwr, Cmid, Cqtr, C3qt, Cmin;
   double Fupr, Flwr, Fmid, Fqtr, F3qt, Fmin;
   double Fold, Cold;
   bool bCheckRepeat;
   int j;

   bCheckRepeat = true;

   //FILE * pFile;
   //pFile = fopen("McCammon.out", "w");   

   Cobs = m_pC[i];
   qobs = m_pq[i];
   wc = m_pWc[i];
   wq = m_pWq[i];

   /*------------------------------------------
   Assign and evaluate initial points, these
   subdivide the domain into 4 quadrants.
   ----------------------------------------- */   
   Cupr = m_Cupr;
   Fupr = F(Cupr, Cobs, qobs, wc, wq);

   Clwr = m_Clwr;
   Flwr = F(Clwr, Cobs, qobs, wc, wq);
  
repeat:
   Cqtr = Clwr + 0.25*(Cupr - Clwr);
   Fqtr = F(Cqtr, Cobs, qobs, wc, wq);

   Cmid = Clwr + 0.50*(Cupr - Clwr);
   Fmid = F(Cmid, Cobs, qobs, wc, wq);

   C3qt = Clwr + 0.75*(Cupr - Clwr);
   F3qt = F(C3qt, Cobs, qobs, wc, wq);

   /* ---------------------------------------
   Debug statements
   fprintf(pFile, "pnt F               C\n");
   fprintf(pFile, "lwr %E  %E\n", Flwr, Clwr);
   fprintf(pFile, "qtr %E  %E\n", Fqtr, Cqtr);
   fprintf(pFile, "mid %E  %E\n", Fmid, Cmid);
   fprintf(pFile, "3qt %E  %E\n", F3qt, C3qt);
   fprintf(pFile, "upr %E  %E\n", Fupr, Cupr);
   fprintf(pFile, "loc it  F_minimum       C_minimum      dC\n");
   ---------------------------------------- */
   
   /* ---------------------------------------------
   Perform bisections. Each bisection will reduce
   the search space by 50%.
   --------------------------------------------- */
   for(j = 0; j < m_MaxIters; j++)
   {
      //is mid-point current minimum?
      if((Fmid <= Fupr) && (Fmid <= Flwr) && (Fmid <= Fqtr) && (Fmid <= F3qt))
      {
         //store new minimum
         Cmin = Cmid;
         Fmin = Fmid;
         //fprintf(pFile, "mid ");

         //shrink domain
         Clwr = Cqtr;
         Flwr = Fqtr;
         Cupr = C3qt;
         Fupr = F3qt;
      }/* end if() */
      //is quarter-point current minimum?
      else if((Fqtr <= Fupr) && (Fqtr <= Flwr) && (Fqtr <= Fmid) && (Fqtr <= F3qt))
      {
         //store new minimum
         Cmin = Cqtr;
         Fmin = Fqtr;
         //fprintf(pFile, "qtr ");

         //shrink domain
         Cupr = Cmid;
         Fupr = Fmid;
         Cmid = Cqtr;
         Fmid = Fqtr;         
      }/* end if() */
      //is three-quarter-point current minimum?
      else if((F3qt <= Fupr) && (F3qt <= Flwr) && (F3qt <= Fmid) && (F3qt <= Fqtr))
      {
         //store new minimum
         Cmin = C3qt;
         Fmin = F3qt;
         //fprintf(pFile, "3qt ");

         //shrink domain
         Clwr = Cmid;
         Flwr = Fmid;
         Cmid = C3qt;
         Fmid = F3qt;
      }/* end if() */
      //is upper-bound current minimum?
      else if((Fupr <= F3qt) && (Fupr <= Flwr) && (Fupr <= Fmid) && (Fupr <= Fqtr))
      {
         //store the new minimum
         Cmin = Cupr;
         Fmin = Fupr;
         //fprintf(pFile, "upr ");

         //shrink domain
         Clwr = C3qt;
         Flwr = F3qt;

         //assign and evaluate new mid-point
         Cmid = Clwr + 0.5*(Cupr - Clwr);
         Fmid = F(Cmid, Cobs, qobs, wc, wq);
      }/* end if() */
      //is lower-bound current minimum?
      else if((Flwr <= F3qt) && (Flwr <= Fupr) && (Flwr <= Fmid) && (Flwr <= Fqtr))
      {
         //store the new minimum
         Cmin = Clwr;
         Fmin = Flwr;
         //fprintf(pFile, "lwr ");

         //shrink domain
         Cupr = Cqtr;
         Fupr = Fqtr;

         //assign and evaluate new mid-point
         Cmid = Clwr + 0.5*(Cupr - Clwr);
         Fmid = F(Cmid, Cobs, qobs, wc, wq);
      }/* end if() */
      //assume mid-point
      else
      {
         //store new minimum
         Cmin = Cmid;
         Fmin = Fmid;
         //fprintf(pFile, "unk ");

         //shrink domain
         Clwr = Cqtr;
         Flwr = Fqtr;
         Cupr = C3qt;
         Fupr = F3qt;
      }/* end if() */

      //assign and evaluate new quarter points
      Cqtr = Clwr + 0.25*(Cupr - Clwr);
      Fqtr = F(Cqtr, Cobs, qobs, wc, wq);

      C3qt = Clwr + 0.75*(Cupr - Clwr);
      F3qt = F(C3qt, Cobs, qobs, wc, wq);

      //fprintf(pFile, "%02d  %E  %E  %E\n", j, Fmin, Cmin, (Cupr - Clwr));
      //fprintf(pFile, "%02d  %E  %E  %E  %E  %E\n", j, Clwr, Cqtr, Cmid, C3qt, Cupr);
   }/* end for() */

   /* -------------------------------------------------
   It's possible that the design space is multi-modal.
   If it is, we may have converged on a sub-optimal 
   point. As a rough guard against such behavior, repeat 
   the bisection search but center the search on Cobs.
   ------------------------------------------------- */
   if(bCheckRepeat == true) //only repeat once
   {      
      bCheckRepeat = false;
      //fprintf(pFile, "Fobs=%E, Fmin=%E, HUGE=%E, NEARLY HUGE=%E\n", Fobs, Fmin, HUGE_VAL, NEARLY_HUGE);
      Clwr = 0.00;
      Flwr = F(Clwr, Cobs, qobs, wc, wq);
      Cupr = 2.00*Cobs;
      Fupr = F(Cupr, Cobs, qobs, wc, wq);
      //save results of initial
      Fold = Fmin;
      Cold = Cmin;                   
      goto repeat;
   }

   if(Fold < Fmin)
   {
      Fmin = Fold; 
      Cmin = Cold;
   }

   // fclose(pFile);

   return Cmin;
}/* end BisectionSearch() */

/******************************************************************************
McCammonSolver::F()

Compute objective function for BiSection Search. Takes five arguments, 
C, observed C, observed q, weight of observed C, and weight of observed q.
******************************************************************************/
double McCammonSolver::F(double C, double Cobs, double qobs, double wc, double wq)
{   
   double f, dq, q;

   q = m_pIso->q(C);
   dq = m_pIso->dqdc(C);

   //equation 10 of McCammon (Am. J. Phys., 1973, v. 5, n. 4, pg. 368)
   f = (Cobs - C)+ (dq*((wq*wq)/(wc*wc))*(qobs - q));

   f *= f;

   //overflow and divide-by-zero can mess up the bisection search
   if(CheckOverflow(f) == true) f = NEARLY_HUGE;

   return f;
}/* end F() */

/******************************************************************************
McCammonSolver::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool McCammonSolver::Initialize(char * pStr)
{
   int i, count;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char * pLine;
   char msg[DEF_STR_SZ];

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginMcCammonMethod") == NULL){ strcat(msg, "BeginMcCammonMethod\n");}
   if(strstr(pStr, "EndMcCammonMethod")   == NULL){ strcat(msg, "EndMcCammonMethod\n");}
   if(strstr(pStr, "BeginConcentrations") == NULL){ strcat(msg, "BeginConcentrations\n");}
   if(strstr(pStr, "EndConcentrations")   == NULL){ strcat(msg, "EndConcentrations\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Orear section
   pTmp = strstr(pStr,"BeginMcCammonMethod");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndMcCammonMethod") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "MaxBisections") != NULL){ sscanf(pLine, "%s %d", pVar, &m_MaxIters);}
   }

   //parse the Concentrations section
   count = 0;
   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndConcentrations") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      count++;
   }
   count--;

   if(count != m_NumOut)
   {
      printf("# of Aqueous/Sorbed Concentrations != # of Total Concentrations\n");
      return false;
   }

   NEW_PRINT("double", m_NumOut);
   m_pWc = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pWq = new double[m_NumOut];

   NEW_PRINT("double", m_NumOut);
   m_pq = new double[m_NumOut];
   MEM_CHECK(m_pq);

   pTmp = strstr(pStr,"BeginConcentrations");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%s %lf %lf %lf %lf", pVar, &(m_pC[i]), &(m_pq[i]), 
                                         &(m_pWc[i]), &(m_pWq[i]));
   }
   return true;
}/* end Initialize() */


