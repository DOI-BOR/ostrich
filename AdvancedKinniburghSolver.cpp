/******************************************************************************
File     : AdvancedKinniburghSolver.cpp
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The AdvancedKinniburghSolver class solves the non-linear Kinniburgh isotherm equation.
See AdvancedKinniburgh.cpp for details.

Version History
03-10-10    lsm   created from KinniburghSolver.h
******************************************************************************/
#include <string.h>

#include "AdvancedKinniburghSolver.h"
#include "Isotherms.h"

#include "Utility.h"
#include "IsoParse.h"
#include "Exception.h"

/******************************************************************************
CTOR

Constructs a AdvancedKinniburghSolver class using the given Isotherm.
******************************************************************************/
AdvancedKinniburghSolver::AdvancedKinniburghSolver(IsothermABC * pIso, double X)
{
   int i;
   m_X = X;
   m_pA = NULL;
   m_pB = NULL;
   m_pD = NULL;
   m_MaxIters = 50;
   
   m_pIso = pIso;

   //share certain data members array with the Isotherm class
   m_pC = pIso->GetPtrToC(&m_NumOut);
   m_pOutFile = pIso->GetPtrToOutFile();

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
      //printf("Kinniburgh() Overflow condition detected and corrected");
   }
   else
   {
      m_Clwr = 0.00;
   }
   m_Cupr *= 2.00;

   IncCtorCount();
} /* end AdvancedKinniburghSolver::CTOR */

/******************************************************************************
DTOR

Destroys a AdvancedKinniburghSolver class.
******************************************************************************/
void AdvancedKinniburghSolver::Destroy(void)
{   
   delete [] m_pA;
   delete [] m_pB;
   delete [] m_pD;
   IncDtorCount();
} /* end AdvancedKinniburghSolver::DTOR */

/******************************************************************************
AdvancedKinniburghSolver::Compute()

Compute output vales and write them to the output file.
******************************************************************************/
void AdvancedKinniburghSolver::Compute(void)
{
   FILE * pFile;
   int i;

   /* -----------------------------------------------------------------
   For each data point, Calculate C that minimizes the following 
   nonlinear equation:
      minimize |C - [Ct/(1+X) - [(S/V)/(1+X)] * q(C)]|
   Since there is  only one equation for each data point, we use a 
   simple "brute-force" bisection method.
   ----------------------------------------------------------------- */
   // perform bisection search
   for(i = 0; i < m_NumOut; i++)
   {
      m_pC[i] = BisectionSearch(i);
   }/* end for() */

   //Utilize isotherm to write out optimal C, q
   m_pIso->Compute();

   //write out Kinniburgh settings (Ct, S, V, and X)
   pFile = fopen(m_pOutFile, "a");
   fprintf(pFile, "\nSolutionMethod AdvancedKinniburgh\n");
   fprintf(pFile, "Max Bisections %d\n", m_MaxIters);
   fprintf(pFile, "Kinniburgh Loss Term %E\n", m_X);

   fprintf(pFile, "\nExperimental Constants\n");
   fprintf(pFile, "i     A(user-defined)  B(user-defined)  D(user-defined)\n");
   for(i = 0; i < m_NumOut; i++)
   {
      fprintf(pFile,"%02d  %E    %E    %E\n",i,m_pA[i], m_pB[i], m_pD[i]);
   }   
   fclose(pFile);
}/* end Compute() */

/******************************************************************************
AdvancedKinniburghSolver::BisectionSearch()

Simple bisection search. Evaluates at least two new points each iteration and 
reduces search space by 50%.

Returns optimial concentration for the i-th observation point.
******************************************************************************/
double AdvancedKinniburghSolver::BisectionSearch(int i)
{
   double Cupr, Clwr, Cmid, Cqtr, C3qt, A, Cmin;
   double Fupr, Flwr, Fmid, Fqtr, F3qt, Fmin, BD;
   int j;
   double Fold, Cold;
   bool bCheckRepeat;

   bCheckRepeat = true;

   //FILE * pFile;
   //pFile = fopen("Kinniburgh.out", "w");   

   A = m_pA[i];
   BD = (m_pB[i]/m_pD[i]);

   /*------------------------------------------
   Assign and evaluate initial points, these
   subdivide the domain into 4 quadrants.
   ----------------------------------------- */   
   Cupr = m_Cupr;
   Fupr = F(Cupr, A, BD);

   Clwr = m_Clwr;
   Flwr = F(Clwr, A, BD);

repeat:
   Cqtr = Clwr + 0.25*(Cupr - Clwr);
   Fqtr = F(Cqtr, A, BD);

   Cmid = Clwr + 0.50*(Cupr - Clwr);
   Fmid = F(Cmid, A, BD);

   C3qt = Clwr + 0.75*(Cupr - Clwr);
   F3qt = F(C3qt, A, BD);

   /* ---------------------------------------
   Debug statements
   fprintf(pFile, "pnt F               C\n");
   fprintf(pFile, "lwr %E  %E\n", Flwr, Clwr);
   fprintf(pFile, "qtr %E  %E\n", Fqtr, Cqtr);
   fprintf(pFile, "mid %E  %E\n", Fmid, Cmid);
   fprintf(pFile, "3qt %E  %E\n", F3qt, C3qt);
   fprintf(pFile, "upr %E  %E\n", Fupr, Cupr);
   fprintf(pFile, "tot %E  %E\n", Ctot, Ctot);
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
         Fmid = F(Cmid, A, BD);
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
         Fmid = F(Cmid, A, BD);
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
      Fqtr = F(Cqtr, A, BD);

      C3qt = Clwr + 0.75*(Cupr - Clwr);
      F3qt = F(C3qt, A, BD);

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
      Flwr = F(Clwr, A, BD);
      Cupr = 2.00*m_pC[i];
      Fupr = F(Cupr, A, BD);
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

   //fclose(pFile);

   return Cmin;  
}/* end BisectionSearch() */

/******************************************************************************
AdvancedKinniburghSolver::F()

Compute objective function for BiSection Search. Takes two arguments, 
C and Ctot.
******************************************************************************/
double AdvancedKinniburghSolver::F(double C, double A, double BD)
{   
   double f;

   /* -----------------------------------------------------------------
   Eqn. 7 of Kinniburgh (Env. Sci. Tech, 1986, v. 20, no. 9, pg. 899)
   It has been generalized to the following form, to allow for 
   flexibility in the definition of experimental constants:
      C = A/(1+X) - [(B/D)/(1+X)]q(C) --> f = C - A/(1+X) + [(B/D)/(1+X)]q(c)
   ----------------------------------------------------------------- */
   f = C - (A/(1.00+m_X)) + ((BD/(1.00+m_X)) * m_pIso->q(C));
   f *= f;

   //overflow and divide-by-zero can mess up the bisection search
   if(CheckOverflow(f) == true) f = NEARLY_HUGE;

   return f;
}/* end F() */

/******************************************************************************
AdvancedKinniburghSolver::Initialize()

Initialize parameters and output arrays, using input file.

Returns true if initialization succeeds, false otherwise.
******************************************************************************/
bool AdvancedKinniburghSolver::Initialize(char * pStr)
{
   int i, count;
   char * pTmp;
   char pVar[DEF_STR_SZ];
   char * pLine;
   char msg[DEF_STR_SZ];

   //check tokens
   msg[0] = (char)NULL;
   if(strstr(pStr, "BeginAdvancedKinniburghMethod") == NULL){ strcat(msg, "BeginAdvancedKinniburghMethod\n");}
   if(strstr(pStr, "EndAdvancedKinniburghMethod")   == NULL){ strcat(msg, "EndAdvancedKinniburghMethod\n");}
   if(strstr(pStr, "BeginExperimentalConstants") == NULL){ strcat(msg, "BeginExperimentalConstants\n");}
   if(strstr(pStr, "EndExperimentalConstants")   == NULL){ strcat(msg, "EndExperimentalConstants\n");}
   if(msg[0] != (char)NULL)
   {
      printf("The following tokens are missing:\n%s", msg);
      return false;
   }

   //parse the Kinniburgh section
   pTmp = strstr(pStr,"BeginAdvancedKinniburghMethod");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndAdvancedKinniburghMethod") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      if(strstr(pLine, "MaxBisections") != NULL){ sscanf(pLine, "%s %d", pVar, &m_MaxIters);}
   }

   //parse the Total Concentrations section
   count = 0;
   pTmp = strstr(pStr,"BeginExperimentalConstants");
   pTmp = ISO_GetLine(pTmp, &pLine);
   while(strcmp(pLine, "EndExperimentalConstants") != 0)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      count++;
   }
   count--;

   if(count != m_NumOut)
   {
      printf("# of Aqueous/Sorbed Concentrations != # of Experimental Constants\n");
      return false;
   }

   NEW_PRINT("double", m_NumOut);
   m_pA = new double[m_NumOut];
   
   NEW_PRINT("double", m_NumOut);
   m_pB = new double[m_NumOut];
   
   NEW_PRINT("double", m_NumOut);
   m_pD = new double[m_NumOut];
   MEM_CHECK(m_pD);

   pTmp = strstr(pStr,"BeginExperimentalConstants");
   pTmp = ISO_GetLine(pTmp, &pLine);
   for(i = 0; i < m_NumOut; i++)
   {
      pTmp = ISO_GetLine(pTmp, &pLine);
      sscanf(pLine, "%lf %lf %lf", &(m_pA[i]), &(m_pB[i]), &(m_pD[i]));
   }
   return true;
}/* end Initialize() */

