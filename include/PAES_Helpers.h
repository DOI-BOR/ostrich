/******************************************************************************
File       : PAES_Helpers.h
Author     : L. Shawn Matott (adpated from Antonio Jesus Nebro Urbaneja)
Copyrights : 2017, L. Shawn Matott (adaptation into OSTRICH)
           : 2004, Antonio Jesus Nebro Urbaneja (original implementation)
           : 2003, Carlos A. Coello Coello
           : 2003, Francisco Luna Valero

This file includes some helper functions and classes for PAES

Version History
12-29-17    lsm   created file
******************************************************************************/
#ifndef PAES_HELPERS_H
#define PAES_HELPERS_H

#include <math.h>
#include "MyTypes.h"

#ifdef __OLD_CPP__
   #include <fstream.h>
   #include <string.h>
   #include <limits.h>
   #include <float.h>
   typedef int bool;
   const bool true  = 1;
   const bool false = 0;
#else
   #include <string>
   #include <fstream>
   #include <climits>
   #include <cfloat>

   const long PAES_MAX_INT = LONG_MAX;
   const long PAES_MIN_INT = LONG_MIN;

   const double PAES_MAX_REAL = NEARLY_HUGE;
   const double PAES_MIN_REAL = -PAES_MAX_REAL;
#endif

#include <stdlib.h>
#include <stdio.h>

enum PAES_MutationOperator {PAES_BIT_FLIP, PAES_RANDOM, PAES_POLYNOMIAL, PAES_UNIFORM };
enum PAES_VariableType {PAES_BINARY, PAES_BINARY_REAL, PAES_REAL, PAES_INTEGER};

/******************************************************************************
class PAES_Random
******************************************************************************/
class PAES_Random 
{
   public:
      int m_Seed;
      int m_isinit;

      float  m_Rseed; //!< Random numbers seed
      double m_oldrand[55]; //!< Array of 55 random numbers
      int    m_jrand; //!< Current random number
      double m_rndx2; //!< Variable used with random normal deviate
      int    m_rndcalcflag; //!< Variable used with random normal deviate

      PAES_Random(void);
      ~PAES_Random(void);

      //functions from DrC.h and DrC.cpp
      void warmup_random(float random_seed);
      float rndreal(float lo, float hi);
      int rnd(int low, int high);
      float randomperc(void);
      double randomnormaldeviate(void);
      void randomize(void);
      double noise(double mu, double sigma);
      void initrandomnormaldeviate(void);
      int flip(float prob);
      void advance_random(void);

      //functions from RandEvolutiva.h and RandEvolutiva.cpp
      double randreal2(void);
      double Gauss(double sigma);
      double N(double m, double sigma);
      void initrandom(int seed);
      int randint(int lo, int hi);
}; /* end class PAES_Random */

#endif /* PAES_HELPERS_H */
