/******************************************************************************
File      : FortranSupportUtilities.h
Author    : L. Shawn Matott
Copyright : 2011, L. Shawn Matott

C-style routines to help with various pre- and post-processing tasks.

Version History
03-17-11    lsm   created
******************************************************************************/

#ifndef FORTRAN_SUPPORT_UTILITIES
#define FORTRAN_SUPPORT_UTILITIES
#include "MyHeaderInc.h"
extern "C" {
int TestSupportUtilities(void);
bool ExtractParameter(const char * name, const char * tpl, const char * inp, bool bFixed, double * val, char ** name_list, int num_params);
bool WriteFixedFormat(FILE * pFile, double val, const char * fmt);
bool GetFixedFormatValAsStr(char * valStr, double val, char * fmt);
int GetMaxLineSize(FILE * pFile);
}
#endif //FORTRAN_SUPPORT_UTILITIES

