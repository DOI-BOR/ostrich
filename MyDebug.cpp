/******************************************************************************
File     : MyDebug.cpp
Author   : L. Shawn Matott
Copyright: 2014, L. Shawn Matott

Macros to aid in debugging.

Version History
10-16-2014 lsm   created from MyDebug.h
******************************************************************************/
#include <stdio.h>

#include "MyDebug.h"

void DBG_PRINT(const char *a)
{
#if DBG_LEVEL >= 1
   printf("%s\n", (char *)a);
#endif
}

void NEW_PRINT(const char *a, int b)
{
#if DBG_LEVEL >= 2
   printf("new %s[%d]\n", (char *)a, b);
#endif
}

