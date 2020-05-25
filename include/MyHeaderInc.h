/******************************************************************************
File     : MyHeaderInc.h
Author   : L. Shawn Matott
Copyright: 2016, L. Shawn Matott

A universal include file for header files. To help avoid cyclic declarations 
this is the only header file (apart from parent abstract base class definitions) 
that should be included by another header file.

Version History
11-30-2016  lsm   created
******************************************************************************/
#ifndef MY_HEADER_INC_H
#define MY_HEADER_INC_H

#include <stdio.h>    // various functions use FILE * and other stdio types
#include <stdlib.h>   // atof() and atoi() are inlined in many class decs
#include <math.h>     // math functions are inlined in many class declarations
#include "MyTypes.h"  // various constants, typedefs, and data structures
#include "MyDebug.h"  // debug print statements --- inlined in many class functions
#include "GeometryUtility.h" // Point, Circle, etc. data structs

/* macros for access and chdir */
#ifdef WIN32
  #include <direct.h>
  #include <io.h>
  #define MY_ACCESS _access
  #define MY_CHDIR _chdir
#else
  #include <unistd.h>
  #define MY_ACCESS access
  #define MY_CHDIR chdir
#endif

// defined in Exception.cpp but inlined in many CTORS and DTORS
void IncCtorCount(void);
void IncDtorCount(void);

#endif /* MY_HEADER_INC_H */



