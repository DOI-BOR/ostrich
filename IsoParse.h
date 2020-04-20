/******************************************************************************
File     : IsoParse.h
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

Parsing routines.

Version History
03-16-04    lsm   Created
03-24-04    lsm   Added functions to support IsoFit program, 
                  which drives Ostrich and internal Isotherm() model
08-31-05    lsm   Added support for Ranges and Swarm Parameters sections
01-01-07    lsm   Added step field to the IsoParamList so that each parameter
                  can have it's own user-specified FD step-size.
******************************************************************************/
#ifndef ISO_PARSE_H
#define ISO_PARSE_H

#include "MyHeaderInc.h"

//forward decs
class ParameterGroup;
class ObservationGroup;

#define ISO_TPL_FILE  "IsothermIn.tpl"
#define ISO_IN_FILE   "IsothermIn.txt"
#define ISO_OSTIN_FILE "ostIn.txt"
#define ISO_OSTOUT_FILE "OstOutput0.txt"
#define ISO_OSTMDL_FILE "OstModel0.txt"
#define ISO_PSOOUT_FILE "OstOutputPSO"
#define ISO_PSOMDL_FILE "OstModelPSO"
#define ISO_OUT_FILE  "IsothermOut.txt"

extern "C" {

#ifdef ISOFIT_BUILD
   int Ostrich(int argc, StringType argv[]);
#endif

typedef enum ISOFIT_SOLVER_ENUM
{
   ISOTHERM_METHOD    = 0,
   OREAR_METHOD       = 1,
   MCCAMMON_METHOD    = 2,
   KINNIBURGH_METHOD  = 3,
   TOTAL_ERROR_METHOD = 4,
   ADV_KINNIBURGH_METHOD =5
}IsoFitSolverType;

typedef struct ISO_PARAM_LIST
{
   char name[20], txin[10], txost[10], txout[10];
   double init, upr, lwr, step;
   struct ISO_PARAM_LIST * pNxt;
}IsoParamList;

//a structure that stores all kinds of data that is read from the
//IsoFit input file and passed around to various subroutines
typedef struct ISO_GLOB_STRUCT
{
   bool bLumpedQ0;
   bool bFitAll;
   bool bHoldObs;
   bool bHoldParams;
   int numObs;
   int debug;
   int popSize;
   int maxGens;
   int maxBisections;
   IsoFitSolverType method;
   char isoStr[DEF_STR_SZ];
   char solStr[DEF_STR_SZ];
   char lumpStr[DEF_STR_SZ];
   double * conc;
   double * sorb;
   double * wsorb; //sorbed concentration weights
   double * wconc; //aqueous concentration weights
   double * expA; //experimental constants
   double * expB; //experimental constants
   double * expD; //experimental constants
}IsoGlobStruct;

int  ISO_GetFileSize(const char * file);
void ISO_FileToStr(const char * file, char * pStr, int size);
char * ISO_GetLine(char * pStr, char ** pLine);

void ISO_ReadIsoFitFile(IsoGlobStruct * pArgs);

void ISO_CreateTemplateFile(IsoGlobStruct * pArgs);

void ISO_CreateOstrichFile(const char * pProgType, IsoParamList * pList, bool stats, 
                           IsoGlobStruct * pArgs);

void ISO_CreateParamList(IsoParamList * pList, IsoGlobStruct * pArgs);
void ISO_RefreshParamList(IsoParamList * pList);
void ISO_DestroyIsoParamList(IsoParamList * pList);
int Isotherm(bool bSave);
int Orear(void);
int McCammon(bool bSave);
int Kinniburgh(bool bSave);
int AdvancedKinniburgh(void);
char * ISO_GetRangesSection(void);
void ISO_GetSwarmParams(char * pStr, IsoGlobStruct * pArgs);
IsoFitSolverType ISO_GetMethod(char * pStr);
void ISO_GetSolutionSettings(char * pStr, IsoGlobStruct * pArgs);
}

void DisklessIsotherm(ParameterGroup * pgroup, ObservationGroup * ogroup);
void DisklessMcCammon(ParameterGroup * pgroup, ObservationGroup * ogroup);
void DisklessKinniburgh(ParameterGroup * pgroup, ObservationGroup * ogroup);

#endif /* ISO_PARSE_H */
