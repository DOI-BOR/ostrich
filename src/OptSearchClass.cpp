/******************************************************************************
File     : OptSearchClass.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

The OptSearchClass is used to perform one-dimensional searches for the optimum
step size of a given search direction.

Version History
08-25-03    lsm   created file
08-18-04    lsm   Added metrics collection and reporting, memory fragmentation
                  fixes, Golden Section search modifications
01-01-07    lsm   OptSearchClass now uses abstract model base class (ModelABC).
******************************************************************************/
#include <math.h>
#include <string.h>

#include "ModelABC.h"
#include "OptSearchClass.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"

const double gCGOLD=(0.5 * (3.00 - sqrt(5.00)));
const double gGOLD=(2.0 - (0.5 * (3.00 - sqrt(5.00))));
const int MaxBrentIters = 36;

/******************************************************************************
CTOR

Initializes everything using m_pModel input.
******************************************************************************/
OptSearchClass::OptSearchClass(ModelABC * pModel)
{
   int i;

   m_pModel = pModel;
   m_NumParams = pModel->GetParamGroupPtr()->GetNumParams();

   NEW_PRINT("double", m_NumParams);
   m_pDir = new double[m_NumParams];
   MEM_CHECK(m_pDir);

   NEW_PRINT("double", m_NumParams);
   m_pStepPoint = new double[m_NumParams];
   MEM_CHECK(m_pStepPoint);

   NEW_PRINT("double", m_NumParams);
   m_pAlphaPoint = new double[m_NumParams];
   MEM_CHECK(m_pAlphaPoint);

   NEW_PRINT("double", m_NumParams);
   m_pStartPoint = new double[m_NumParams];
   MEM_CHECK(m_pStartPoint);

   NEW_PRINT("MinBracketStruct", 1);
   m_pMinBrack = new MinBracketStruct;
   MEM_CHECK(m_pMinBrack);

   NEW_PRINT("MinPtStruct", 1);
   m_pMinPt = new MinPtStruct;
   MEM_CHECK(m_pMinPt);

   m_Step = 0.00;
   m_SearchConvVal = 1E-4;
   m_BoundMinCount = 0;
   m_GoldSectCount = 0;
   m_BrentCount    = 0;
   m_SearchType    = GSECT_SEARCH;
   
   //init. direction
   for(i = 0; i < m_NumParams; i++)
   {
      m_pDir[i] = 0.00;
   }/* end for() */

   InitFromFile(GetInFileName());

   IncCtorCount();
}/* end CTOR */

/******************************************************************************
Destroy()

Frees up memory.
******************************************************************************/
void OptSearchClass::Destroy(void)
{
   delete [] m_pDir;
   delete [] m_pStepPoint;
   delete [] m_pAlphaPoint;
   delete [] m_pStartPoint;
   delete m_pMinBrack;
   delete m_pMinPt;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WriteMetrics()

Reports on the setup of the Search Class and also various run-time metrics.
******************************************************************************/
void OptSearchClass::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nOne Dimensional Search Metrics\n");
   fprintf(pFile, "Search Type          : ");
   if(m_SearchType == GSECT_SEARCH){ fprintf(pFile, "Golden Section\n");}   
   else{ fprintf(pFile, "Brent\n");}
   fprintf(pFile, "Convergence Val      : %E\n", m_SearchConvVal);
   fprintf(pFile, "Bound Min Evals      : %d\n", m_BoundMinCount);
   fprintf(pFile, "Golden Section Evals : %d\n", m_GoldSectCount);
   fprintf(pFile, "Brent Evals          : %d\n", m_BrentCount);
}/* end WriteMetrics() */

/******************************************************************************
InitFromFile()

Reads configuration parameters from input file.
******************************************************************************/
void OptSearchClass::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];

   //read in configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open 1-D search config. file. Using Defaults");
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "Begin1dSearch", pFileName) == true)
   {
      FindToken(pFile, "End1dSearch", pFileName);
      rewind(pFile);

      FindToken(pFile, "Begin1dSearch", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "End1dSearch") == NULL)
      {
         if(strstr(line, "1dSearchConvergeVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_SearchConvVal);
         }/*end if() */         
         else if(strstr(line, "1dSearchMethod") != NULL)
         {
            MyStrLwr(line);
            if(strstr(line, "brent") != NULL)
            {
               m_SearchType = BRENT_SEARCH;
            }
            else if(strstr(line, "goldensection") != NULL)
            {
               m_SearchType = GSECT_SEARCH;
            }
            else
            {
               strcpy(tmp, "Unknown Search Method");
               LogError(ERR_FILE_IO, tmp);
            }
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
}/* end InitFromFile() */

/******************************************************************************
CalcF()

Compute objective function, using a step size of <alpha>.

NOTE: After completion of the calculation, the model is NOT rerun 
at the initial location, therefore this routine will leave the system in an 
inconsistent state. Routines that call CalcF() must handle the restoration of 
the system.

If a better minimum than fmin is found, fmin and xmin are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
double OptSearchClass::CalcF(double alpha, double * fmin, double * xmin)
{
   ParameterGroup * pGroup;
   double upr, lwr, Finit;
   int i;
   double F, old_pi;

   //backup initial location of design point
   pGroup = m_pModel->GetParamGroupPtr();
   pGroup->ReadParams(m_pStartPoint);
   Finit = m_pModel->GetObjFuncVal();
   pGroup->ReadParams(m_pAlphaPoint);
   
   /*---------------------------------------------------------
   Adjust each design parameter by step size, taking care to
   avoid stepping out of the side constraints of any given
   parameter
   ---------------------------------------------------------*/
   for(i = 0; i < m_NumParams; i++)
   {
      old_pi = m_pAlphaPoint[i];

      m_pAlphaPoint[i] = m_pAlphaPoint[i] + alpha * m_pDir[i];

      //if it's out of bounds, move half the distance to upr/lwr
      upr = pGroup->GetParamPtr(i)->GetUprBnd();
      lwr = pGroup->GetParamPtr(i)->GetLwrBnd();
      if(m_pAlphaPoint[i] > upr)
	  {
		  m_pAlphaPoint[i] = (upr+old_pi)/2.00;
	  }
      if(m_pAlphaPoint[i] < lwr)
	  {
		  m_pAlphaPoint[i] = (old_pi+lwr)/2.00;
	  }
   }
   //run model at new location
   pGroup->WriteParams(m_pAlphaPoint);
   F = m_pModel->Execute();
   /* update optimal, if appropriate */
   if(F < *fmin)
   {
	   *fmin = F;
	   pGroup->ReadParams(xmin);
   }

   //semi-restore model (for next time around)
   pGroup->WriteParams(m_pStartPoint);
   m_pModel->SetObjFuncVal(Finit);

   return F;
}/* end CalcF() */

/******************************************************************************
CalcStepSize()

Compute the optimum step size in the direction given by <pDir>.

If a better minimum than fmin is found, fmin and xmin are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
double OptSearchClass::CalcStepSize(Unchangeable1DArray pDir, double * fmin, double * xmin)
{   
   int i;
   MinBracketStruct * pMbrak;
   MinPtStruct * pMin;
   double Finit, Fcur;

   //store current setting
   m_pModel->GetParamGroupPtr()->ReadParams(m_pStepPoint);
   Finit = m_pModel->GetObjFuncVal();

   //store direction
   for(i = 0; i < m_NumParams; i++)
   {
      m_pDir[i] = pDir[i];
   }/* end for() */
   
   //bound minimum
   pMbrak = BracketMinimum(-1.00, 1.00, fmin, xmin);

   /* Golden Section Method */
   if(m_SearchType == GSECT_SEARCH){ pMin = GoldSect(pMbrak, fmin, xmin);}
   /* Brents method */
   else{ pMin = Brent(pMbrak, fmin, xmin);}

   m_Step = pMin->x;   

   //restore initial setting, ensuring that model is consistent
   m_pModel->GetParamGroupPtr()->WriteParams(m_pStepPoint);
   Fcur = m_pModel->Execute();
   m_BoundMinCount++;
   //check model consistency
   if(Fcur != Finit)
   {
      LogError(ERR_MODL_EXE, "CalcStepSize() caused model to be inconsistent");
   }
   
   return (m_Step);
}/* CalcStepSize() */

/******************************************************************************
BracketMinimum()

Brackets the minimum using a standard Golden Section algorithm. 

Returns: a MinBracket structure containing points a, b and c along with f(a),
f(b), and f(c).

If a better minimum than fmin is found, fmin and xmin are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
MinBracketStruct * OptSearchClass::BracketMinimum(double a, double b, double * fmin, double * xmin)
{
   double ulim, u, r, q, fu, dum, Fcur, denom;
   double fa, fb, c, fc;

   //initial values
   Fcur = m_pModel->GetObjFuncVal();

   a = LimitStepSize(a, m_pDir);
   fa = CalcF(a, fmin, xmin);
   m_BoundMinCount++;   

   b = LimitStepSize(b, m_pDir);
   fb = CalcF(b, fmin, xmin);
   m_BoundMinCount++;  
   
   if(fb >  fa) //switch roles of a and b
   {
     dum = a;
     a = b;
     b = dum;

     dum = fa;
     fa = fb;
     fb = dum;
   }/* end if() */

   //check to see if initial fb and fa already bracket a minimum
   if((fb >= Fcur) && (fa >= Fcur))
   {
      c = b;
      fc = fb;
      b = 0.00;
      fb = Fcur;
   }/* end if() */
   else
   {
      c = b + (gGOLD * (b - a));
	  c = LimitStepSize(c, m_pDir);
      fc = CalcF(c, fmin, xmin);
      m_BoundMinCount++;

      while(fb > fc)
      {
         r = (b - a)*(fb - fc);
         q = (b - c)*(fb - fa);
         denom = MyMax(fabs(q-r),NEARLY_ZERO);
         if((q-r) < 0) denom = -denom;
         u = (b) - ((b - c)*q - (b - a)*r)/(2.0*denom);
         ulim = b + 100*(c - b);
         //test various possibilities
         if((b - u)*(u - c) > 0.0)
         {
			u = LimitStepSize(u, m_pDir);
            fu = CalcF(u, fmin, xmin);
            m_BoundMinCount++;

            if(fu <= fc) //min. between b and c
            {
               a = b;
               b = u;
               fa = fb;
               fb = fu;
               break;
            }/* end if() */
            else if(fu >= fb) //min. between a and u
            {
               c = u;
               fc = fu;
               break;
            } /* end else if() */
            u = (c) + gGOLD*(c - b);
			u = LimitStepSize(u, m_pDir);
            fu = CalcF(u, fmin, xmin);
            m_BoundMinCount++;
         }/* end if() */
         else if((c - u)*(u - ulim) > 0.0)
         {
			u = LimitStepSize(u, m_pDir);
            fu = CalcF(u, fmin, xmin);
            m_BoundMinCount++;
            if(fu < fc)
            {
               dum = c + gGOLD*(c-b);
               b = c;
               c = u;
               u = dum;

               dum = CalcF(u, fmin, xmin);
               fb = fc;
               fc = fu;
               fu = dum;

               m_BoundMinCount++;
            }/* end if() */
         }/* end else if() */
         else if((u - ulim)*(ulim - c) >= 0.0)
         {
            u = ulim;
			u = LimitStepSize(u, m_pDir);
            fu = CalcF(u, fmin, xmin);
            m_BoundMinCount++;
         }/* end else if() */
         else
         {
            u = (c) + gGOLD*(c - b);
			u = LimitStepSize(u, m_pDir);
            fu = CalcF(u, fmin, xmin);
            m_BoundMinCount++;
         }/* end else() */         
         a = b;
         b = c;
         c = u;

         fa = fb;
         fb = fc;
         fc = fu;
      }/* end while() */
   } /* end else() */

   m_pMinBrack->a = a;
   m_pMinBrack->b = b;
   m_pMinBrack->c = c;
   m_pMinBrack->fa = fa;
   m_pMinBrack->fb = fb;
   m_pMinBrack->fc = fc;
      
   return m_pMinBrack;
}/* end BracketMinimum() */

/******************************************************************************
GoldSect()

Calculate the minimum using the golden section method, as described by 
Vanderplaats (page 57-60).

Returns location of the min. point.

If a better minimum than minf is found, minf and minp are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
MinPtStruct * OptSearchClass::GoldSect(MinBracketStruct * pBrack, double * minf, double * minp)
{   
   double x1, x2, F1, F2, xL, xU, Fmin, xmin;
   int i, its;

   //determine initial min
   Fmin = pBrack->fa; xmin = pBrack->a;
   if(Fmin > pBrack->fb){ Fmin = pBrack->fb; xmin = pBrack->b;}
   if(Fmin > pBrack->fc){ Fmin = pBrack->fc; xmin = pBrack->c;}

   xL = pBrack->a;   
   xU = pBrack->c;   

   /*---------------------------------------------------------------------
   Calculate the number of iterations required to reduce the initial width
   to the desired convergence value.
   ---------------------------------------------------------------------*/
   its = (int)((log10(m_SearchConvVal/fabs(xL-xU))/log10(1.00-gCGOLD)) + 1.00);

   Write1dSearch(WRITE_GSECT, its+2);
   x1 = ((1.00 - gCGOLD) * xL) + (gCGOLD * xU);
   x1 = LimitStepSize(x1, m_pDir);
   F1 = CalcF(x1, minf, minp);
   if(Fmin > F1){ Fmin = F1; xmin = x1;}
   m_GoldSectCount++;
   Write1dSearch(1, 0);

   x2 = ((1.00 - gCGOLD) * xU) + (gCGOLD * xL);
   x2 = LimitStepSize(x2, m_pDir);
   F2 = CalcF(x2, minf, minp);
   if(Fmin > F2){ Fmin = F2; xmin = x2;}
   m_GoldSectCount++;
   Write1dSearch(2, 0);
      
   //Golden Section iteration   
   for(i = 0; i < its; i++)
   {
      Write1dSearch(i+3, 0);
      if(F1 > F2)
      {
         xL = x1;
         x1 = x2;
         F1 = F2;
         x2 = ((1.00 - gCGOLD) * xU) + (gCGOLD * xL);
		 x2 = LimitStepSize(x2, m_pDir);
         F2 = CalcF(x2, minf, minp);
         if(Fmin > F2){ Fmin = F2; xmin = x2;}
         m_GoldSectCount++;
      }/* end if() */
      else /* (F2 > F1) */
      {
         xU = x2;
         x2 = x1;
         F2 = F1;
         x1 = ((1.00 - gCGOLD) * xL) + (gCGOLD * xU);
		 x1 = LimitStepSize(x1, m_pDir);
         F1 = CalcF(x1, minf, minp);
         if(Fmin > F1){ Fmin = F1; xmin = x1;}
         m_GoldSectCount++;      
      }/* end else() */
   }/* end while() */
         
   m_pMinPt->x  = xmin; 
   m_pMinPt->fx = Fmin;

   Write1dSearch(WRITE_ENDED, 0);   
   return m_pMinPt;
}/* end GoldSect() */

/******************************************************************************
Brent()

Finds the minimum, given the bracketing triplet stored in MinBracketStruct arg.
Uses the well-known Brent's method. 

Returns: a MinPtStruct containing the mininum point and the value of the 
objective function at that point.

If a better minimum than minf is found, minf and minp are replaced with the
new minimum value and corresponding parameter values.
******************************************************************************/
MinPtStruct * OptSearchClass::Brent(MinBracketStruct * pBrack, double * minf, double * minp)
{
   int iter;
   double ax, bx, cx, tol, xmin, fmin;
   double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;

   ax = pBrack->a;
   bx = pBrack->b;
   cx = pBrack->c;
   tol = m_SearchConvVal;

   if(ax < cx){ a = ax;} else { a = cx;}
   if(ax > cx){ b = ax;} else { b = cx;}
   x=w=v=bx;
   fw=fv=fx=pBrack->fb;   

   Write1dSearch(WRITE_BRENT, MaxBrentIters);
   for(iter = 1; iter <= MaxBrentIters; iter++)
   {
      Write1dSearch(iter, 0);
      xmin = x;
      fmin = fx;

      xm = 0.5*(a + b);
      tol1 = tol*fabs(x) + NEARLY_ZERO;
      tol2 = 2.0*tol1;
      if(fabs(x - xm) <= (tol2 - 0.5*(b - a)))
      {
         //converged on minimum, so break out of for loop...
         break;
      }/* end if() */
      if(fabs(e) > tol1)
      {
         r = (x - w)*(fx - fv);
         q = (x - v)*(fx - fw);
         p = (x - v)*q - (x - w)*r;
         q = 2.0*(q - r);
         if(q > 0.0){p = -p;}
         q = fabs(q);
         etemp = e;
         e = d;
         if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b - x))
         {  
            if(x >= xm){e = (a - x);} else {e = (b - x);} 
            d = gCGOLD * e;            
         }/* end if() */
         else
         {
            d = p/q;
            u = x + d;
            if(u - a < tol2 || b - u < tol2) 
            {
               if((xm - x) < 0) d = -fabs(tol1);
               else d = fabs(tol1);
            }
         }/* end else() */
      }/* end if() */      
      else
      {  
         if(x >= xm){e = (a - x);} else {e = (b - x);}       
         d = gCGOLD * e;                  
      }
      if(fabs(d) >= tol1)
      { 
         u = x + d;
      } 
      else
      {
         if(d >= 0)
         {
            u = x + tol1;
         }
         else
         {
            u = x - tol1;
         }
      }

	   u = LimitStepSize(u, m_pDir); //prevent step size from exceeding parameter bounds
      fu = CalcF(u, minf, minp);
      m_BrentCount++;
      
      if(fu == fx){ break;} //multi-modal?
      if(fu < fx) //found a better minimum
      {
         if(u >= x){a = x;} else {b = x;}
         v = w;
         w = x;
         x = u;

         fv = fw;
         fw = fx;
         fx = fu;
      }/* end if() */
      else
      {
         if(u < x){a = u;} else{b = u;}
         if(fu <= fw || w == x)
         {
            v = w;
            w = u;
            fv = fw;
            fw = fu;
         }/* end if() */
         else if(fu <= fv || v == x || v == w)
         {
            v = u;
            fv = fu;            
         }/* end else if() */
      }/* end else() */
   }/* end for() */

   //if we exceed the max. iterations or stalled at first iteration, revert to trusty Golden Section method
   if((iter > MaxBrentIters) || (iter == 1))
   { 
      Write1dSearch(WRITE_SWTCH, 0);
      m_pMinPt = GoldSect(pBrack, minf, minp);
   }
   else
   {
      Write1dSearch(WRITE_ENDED, 0);
      m_pMinPt->x  = xmin;
      m_pMinPt->fx = fmin;
   }
   
   return m_pMinPt;
} /* end Brent() */

/******************************************************************************
LimitStepSize()

Prevents step size from exceeding parameter bounds.

Returns: The largest allowable step size given the parameter bounds and the
search directions.
******************************************************************************/
double OptSearchClass::LimitStepSize(double alpha, double * pDir)
{
   int i;
   double upr, lwr, pi, opi;

   /*---------------------------------------------------------
   Check each design parameter to see if it can accommodate
   the given step size.
   ---------------------------------------------------------*/
   for(i = 0; i < m_NumParams; i++)
   {
      opi = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetEstVal();
	  pi = opi + alpha * pDir[i];
      //if it's out of bounds, move half the distance to upr/lwr
      upr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetUprBnd();
      lwr = m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetLwrBnd();
      if(pi > upr)
	  {
         alpha = (upr - opi)/pDir[i];
	  }
      if(pi < lwr)
	  {
         alpha = (lwr - opi)/pDir[i];
	  }
   }
   return alpha;
}/* end LimitStepSize() */
