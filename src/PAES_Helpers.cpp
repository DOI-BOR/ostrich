/******************************************************************************
File       : PAES_Helpers.cpp
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)
           : 2003, Carlos A. Coello Coello
           : 2003, Francisco Luna Valero

This file includes some helper functions and classes for PAES

Version History
12-29-17    lsm   created file
******************************************************************************/
#include "PAES_Helpers.h"

#define PAES_NRAND_SAMPLES 5
#define PAES_Uniform randreal2
#define PAES_Gmaximo 5000

/******************************************************************************
PAES_Random CTOR
******************************************************************************/
PAES_Random::PAES_Random(void) 
{
    m_Seed = 0;
    m_isinit = 0;
} /* end PAES_Random::CTOR */

/******************************************************************************
PAES_Random DTOR
******************************************************************************/
PAES_Random::~PAES_Random(void) 
{
    return;
} /* end PEAS_Random DTOR */

/******************************************************************************
PAES_Random::flip()

Flip a biased coin - true if heads 
    prob: The biased probability
******************************************************************************/
int PAES_Random::flip(float prob)
{
   if(randomperc() <= prob)
   {
	   return(1);
   }
   else
   {
      return(0);
   }
}/* end flip() */

/******************************************************************************
PAES_Random::advance_random()
******************************************************************************/
void PAES_Random::advance_random(void)
{
   int j1;
   double new_random;

   for(j1 = 0; j1 < 24; j1++)
   {
      new_random = m_oldrand[j1] - m_oldrand[j1+31];

      if(new_random < 0.0) 
      {
         new_random = new_random + 1.0;
      }
      m_oldrand[j1] = new_random;
   }/* end for() */

   for(j1 = 24; j1 < 55; j1++)
   {
      new_random = m_oldrand[j1] - m_oldrand[j1-24];

      if(new_random < 0.0) 
      {
         new_random = new_random + 1.0;
      }

      m_oldrand[j1] = new_random;
   }/* end for()  */
}/* end advance_random() */

/******************************************************************************
PAES_Random::initrandomnormaldeviate()

Initialization routine for randomnormaldeviate
******************************************************************************/
void PAES_Random::initrandomnormaldeviate(void)
{
   m_rndcalcflag = 1;
}/* end initrandomnormaldeviate() */

/******************************************************************************
PAES_Random::noise()

Normal noise with specified mean & std dev: mu & sigma 
   mu: The mean of the distribution
   sigma: The standard deviation of the distribution
******************************************************************************/
double PAES_Random::noise(double mu, double sigma)
{
   return((randomnormaldeviate()*sigma) + mu);
}/* end noise() */

/******************************************************************************
randomize()

Initialize a batch of random numbers 
******************************************************************************/
void PAES_Random::randomize(void)
{
   int j1;

   for(j1=0; j1<=54; j1++)
   {
      m_oldrand[j1] = 0.0;
   }

   m_jrand = 0;

   warmup_random(m_Rseed);
}/* end randomize() */

/******************************************************************************
randomnormaldeviate()

Random normal deviate after ACM algorithm 267 / Box-Muller Method 
******************************************************************************/
double PAES_Random::randomnormaldeviate(void)
{
   double t, rndx1;

   if(m_rndcalcflag)
   {
      rndx1 = sqrt(- 2.0*log((double) randomperc()));
      t = 6.2831853072 * (double) randomperc();
      m_rndx2 = sin(t);
      m_rndcalcflag = 0;
      return(rndx1 * cos(t));
   }
   else
   {
      m_rndcalcflag = 1;
      return(m_rndx2);
   }
}/* end randomnormaldeviate() */

/******************************************************************************
randomperc()

Fetch a single random number between 0.0 and 1.0 

Fetch a single random number between 0.0 and 1.0 - Subtractive Method 
See Knuth, D. (1969), v. 2 for details 
name changed from random() to avoid library conflicts on some machines
******************************************************************************/
float PAES_Random::randomperc(void)
{
   m_jrand++;

   if(m_jrand >= 55)
   {
      m_jrand = 1;
      advance_random();
   }
   return((float) m_oldrand[m_jrand]); 
}/* end randomperc() */

/******************************************************************************
rnd()

Pick a random integer between low and high
******************************************************************************/
int PAES_Random::rnd(int low, int high)
{
   int i;
   if(low >= high)
   {
      i = low;
   }
   else
   {
      i =(int) (randomperc() * (high - low + 1)) + low;
      if(i > high)
      { 
         i = high;
      }
   }

   return(i);
}/* end rnd() */

/******************************************************************************
rndreal()

Real random number between specified limits
******************************************************************************/
float PAES_Random::rndreal(float lo, float hi)
{
   return((randomperc() * (hi - lo)) + lo);
}/* end rndreal() */

/******************************************************************************
warmup_random()

Get random off and running
******************************************************************************/
void PAES_Random::warmup_random(float random_seed)
{
   int j1, ii;
   double new_random, prev_random;

   m_oldrand[54] = random_seed;
   new_random = 0.000000001;
   prev_random = random_seed;

   for(j1 = 1; j1 <= 54; j1++)
   {
      ii = (21*j1)%54;
      m_oldrand[ii] = new_random;
      new_random = prev_random-new_random;
      if(new_random<0.0) 
      {
         new_random = new_random + 1.0;
      }
      prev_random = m_oldrand[ii];
   }/* end for() */

   advance_random();
   advance_random();
   advance_random();

   m_jrand = 0;
}/* end warmup_random() */

/******************************************************************************
randreal2()
******************************************************************************/
double PAES_Random::randreal2(void)
{
   double result;

   if (!m_isinit) 
   {
      srand(m_Seed);
      m_isinit = 1;
   }/* end if() */
   result = ((double) rand());
   result /= RAND_MAX;

   return (result);
}/* end randreal2() */

/******************************************************************************
Gauss()

Gauss distribution

sigma: The standard deviation

Returns a Gauss distribution value
******************************************************************************/
double PAES_Random::Gauss(double sigma)
{
   double ret_val;
   static double u, x, y, u0, u1, u2;

//L1:
   u = PAES_Uniform();
   u0 = PAES_Uniform();
   if (u >= 0.919544406) 
   {
      goto L2;
   }
   x = (u0 + u * 0.825339283) * 2.40375766 - 2.11402808;
   goto L10;

L2:
   if (u < 0.965487131) 
   {
      goto L4;
   }

L3:
   u1 = PAES_Uniform();
   y = sqrt(4.46911474 - log(u1) * 2.00);
   u2 = PAES_Uniform();
   if (y * u2 > 2.11402808) 
   {
      goto L3;
   }
   goto L9;

L4:
   if (u < 0.949990709) 
   {
      goto L6;
   }

L5:
   u1 = PAES_Uniform();
   y = u1 * 0.273629336 + 1.84039875;
   u2 = PAES_Uniform();
   if ((exp(y * -0.5 * y) * 0.39894228 - 0.443299126 + y * 0.209694057) < (u2 * 0.0427025816)) 
   {
      goto L5;
   }
   goto L9;

L6:
   if (u < 0.925852334) 
   {
      goto L8;
   }

L7:
   u1 = PAES_Uniform();
   y = u1 * 1.55066917 + 0.289729574;
   u2 = PAES_Uniform();
   if ((exp(y * -0.5 * y) * 0.39894228 - 0.443299126 + y * 0.209694057) < (u2 * 0.0159745227)) 
   {
      goto L7;
   }
   goto L9;

L8:
   u1 = PAES_Uniform();
   y = u1 * 0.289729574;
   u2 = PAES_Uniform();
   if ((exp(y * - 0.5 * y) * 0.39894228 - 0.382544556) < (u2 * 0.0163977244)) 
   {
      goto L8;
   }

L9:
   x = y;
   if (u0 >= 0.5) 
   {
      x = -y;
   }

L10:
   ret_val = sigma * x;

   return ret_val;
} /* Gauss() */

/******************************************************************************
N()
******************************************************************************/
double PAES_Random::N(double m, double sigma)
{
   return m + Gauss(sigma);
}/* end N() */

/******************************************************************************
initrandom()

Initialize the random seed

seed: The random seed
******************************************************************************/
void PAES_Random::initrandom(int seed)
{
	m_Seed = seed;
}/* end initrandom() */


/******************************************************************************
randint()
******************************************************************************/
int PAES_Random::randint(int lo, int hi)
{
   return (lo + (int) randreal2() * (hi - lo + 1));
}/* end randint() */

