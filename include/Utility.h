/******************************************************************************
File      : Utility.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

This file contains a bunch of useful c-style routines, ranging from matrix 
mathematics to string manipulation.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
10-27-03    lsm   Removed all input files except model input file.
11-25-03    lsm   Added string reverse function
08-17-04    lsm   Added string replace and string occurence functions.
                  RAM fragmentation fixes.
11-30-04    lsm   Replaced calls to time() with MyTime()
10-19-05    lsm   Added MyRand() and MY_RAND_MAX
01-01-07    lsm   Added support for temporary input files which store
                  copies of the surrogate sections of the input file.
03-03-07    jrc   Added UniformRandom(), GaussRandom(), and iMax()
******************************************************************************/
#ifndef UTILITY_H
#define UTILITY_H

extern "C" {

#define MY_RAND_MAX (0x7FFFFFFF)
unsigned int MyRand(void);
void MyRandCleanup(void);
double MyGaussRand(double m, double s);
double			 UniformRandom(void);
double			 GaussRandom(void);
double * MultiVarGaussRandom(double * x, double ** covar, int n);
double * OverDispMultiVarGaussRandom(double * x, double ** covar, double s, double w, int n);
unsigned int MyTime(void);
unsigned int GetElapsedTime(void);
double GetElapsedTics(void);
unsigned int GetRandomSeed(void);
void ResetRandomSeed(unsigned int seed);
unsigned int ReadRandomSeed(void);
void RestoreRandomSeed(void);
int SampleWithReplacement(int opFlag, int range);

bool AlgIsMultiObjective(void);

char * GetOstExeOut(void);

bool SynchReceives(void);

bool IsWhitespace(char c);
bool IsNumeric(char c);

int ExtractFileName(IroncladString pLine, UnmoveableString pOut);

int ExtractString(IroncladString pLine, UnmoveableString pOut);
int ExtractColString(IroncladString pLine, UnmoveableString pOut, char tok);
int ValidateExtraction(int j, int cur, int last, const char * pF); 
int CheckExtraction(int j, int cur, int last, const char * pF); 
int GetMaxLineSizeInString(char * pStr);
int GetMaxLineSizeInFile(char * fname);

void FindToken(FILE * pFile, IroncladString token, IroncladString pName);

bool CheckToken(FILE * pFile, IroncladString token, IroncladString pName);

bool CheckOverflow(double num);

char * GetNxtDataLine(FILE * pFile, IroncladString pName);
char * GetCurDataLine(void);
void InitDataLine(IroncladString pName);

void GetPreciseValAsStr(UnmoveableString pOut, double x);

void SortInc(Unmoveable1DArray v, int size);
int CompDbl(UnchangeableVoidPtr arg1,UnchangeableVoidPtr arg2);
void MatMult(Ironclad2DArray m1, Ironclad2DArray m2, Unmoveable2DArray mOut, 
             int row1, int row2, int col2);
void VectMult(Ironclad2DArray m, Ironclad1DArray v, Unmoveable1DArray vOut, 
              int rows, int cols);  
double DotProduct(Ironclad1DArray v1, Ironclad1DArray v2, int size);
bool MatInv(Ironclad2DArray m, Unmoveable2DArray inv, int size);
int CholeskyDecomp(Ironclad2DArray A, Unmoveable2DArray L, Unmoveable2DArray LT, int size);

void MyTrim(UnmoveableString pLine);
void MyStrLwr(UnmoveableString pLine);
void MyStrRev(UnmoveableString pLine);
int  MyStrRep(UnmoveableString pStr, IroncladString pFind, IroncladString pRep);
int  MyStrOccur(UnmoveableString pStr, IroncladString pFind);
void MyStrDiff(UnmoveableString pStr1, UnmoveableString pStr2);
void MyStrProtect(UnmoveableString pStr, IroncladString pName);
void  MyStrUnProtect(UnmoveableString pStr, IroncladString pName);

double MyMin(double a, double b);
double MyMax(double a, double b);
int    iMax(const int a, const int b);

char * MyTempName(UnmoveableString pStr);
UnmoveableString GetOstFileName(void);
UnmoveableString GetExeDirName(void);
UnmoveableString GetInFileName(void);
UnmoveableString GetSrgFileName(void);
UnmoveableString GetDynFileName(IroncladString pTok);

void ExecuteCommandLine(IroncladString cmd, bool isRead, IroncladString fileName, IroncladString paramName);
void ConvertToASCII(void);
ProgramType GetProgramType(void);
void SetProgramType(ProgramType progVal);
ProgramType ReadProgramType(void);
void SetOstExePath(char * pPath);
char * GetOstExePath(void);

int GetNumOutputFiles(char * prefix, char * suffix);
void GetOutputFiles(int num, char ** pList, char * prefix, char * suffix);
int SimpleWarmStart(int np, double * best);
double GetBestObjFunc(int np);
bool IsNonDominated(double * pF, int nObj);

}
#endif /* UTILITY_H */

