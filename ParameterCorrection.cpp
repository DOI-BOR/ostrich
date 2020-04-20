/******************************************************************************
File     : ParameterCorrection.cpp
Author   : L. Shawn Matott
Copyright: 2012, L. Shawn Matott

Intereface for external parameter correction algorithm.

Version History
09-18-12    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ParameterCorrection.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "ResponseVarGroup.h"
#include "RespVarABC.h"
#include "FilePair.h"

#include "Utility.h"
#include "Exception.h"

/******************************************************************************
 default CTOR
******************************************************************************/
ParameterCorrection::ParameterCorrection(ParameterGroup * pGroup)
{
   FilePair * pFilePair;
   char * line;
   char tmp1[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char tmp3[DEF_STR_SZ];   
   FILE * pInFile;
   int i;
   int j;
   bool quoteWrap; //if true, then executable must be wrapped in quotes

   //RegisterCorrectorPtr(this);

   IroncladString inFileName = GetInFileName();
   UnmoveableString pDirName = GetExeDirName();

   m_ExecCmd = NULL;
   m_pCorrections = NULL;
   m_FileList = NULL;
   m_NumCorrections = 0;
   
   pInFile = fopen(inFileName, "r");
   if(pInFile == NULL)
   {
      FileOpenFailure("ParameterCorrection::CTOR", inFileName);
   }/* end if() */

   //check for critical entries, entries which have no reasonable defaults
   FindToken(pInFile, "BeginParameterCorrection", inFileName);
   FindToken(pInFile, "EndParameterCorrection", inFileName);
   rewind(pInFile);
   FindToken(pInFile, "BeginParameterCorrection", inFileName);   

   /*
   ------------------------------------------------------------------
   Read in the corrector executable, modifying so that output is 
   redirected to a file.
   ------------------------------------------------------------------
   */
   do
   {
      line = GetNxtDataLine(pInFile, inFileName);
      if(feof(pInFile)) break;

      /*--------------------------------------------------------
      Read in executable, taking care to preserve full path, 
      even in the presence of long and space-separated filenames.
   -  -------------------------------------------------------*/   
      if(strncmp(line, "Executable", strlen("Executable")) == 0)
      {
         i = ExtractString(line, tmp2);
         i = ValidateExtraction(i, 1, 1, "ParameterCorrection()");
         i = ExtractFileName(&(line[i]), tmp1);

         //must wrap in quotes if there is whitespace in the execuable path
         quoteWrap = false;
         j = (int)strlen(tmp1);
         for(i = 0; i < j; i++){ if(tmp1[i] == ' ') {quoteWrap = true;}}     
         if(quoteWrap == true) { tmp1[j++] = '"';}
         tmp1[j] = (char)NULL;   
         if(quoteWrap == true) 
         { 
            MyStrRev(tmp1);
            tmp1[j++] = '"';
            tmp1[j] = (char)NULL;
            MyStrRev(tmp1);
         }/* end if() */

         //make sure the executable exists
         strcpy(tmp2, tmp1);

         if(tmp2[0] == '"')
         {
            tmp2[0] = ' ';
            tmp2[strlen(tmp2)-1] = ' ';
            MyTrim(tmp2);
         }
         if(MY_ACCESS(tmp2, 0 ) == -1)
         {
            sprintf(tmp1, "Parameter correction executable (|%s|) not found", tmp2);
            LogError(ERR_FILE_IO, tmp1);
            ExitProgram(1);
         }

         #ifdef WIN32 //windows version
            strcat(tmp1, " > OstParameterCorrectionOut.txt");
         #else //Linux (bash, dash, csh)
            strcpy(tmp2, tmp1);
            // '>&' redircts both output and error
            strcat(tmp2, " > OstParameterCorrectionOut.txt 2>&1"); 
            strcpy(tmp1, tmp2);
         #endif
         
         SetExecCmd(tmp1);
      }/* end if() */
      /*
      ---------------------------------------------------------------
      Read in the 'File Pairs': a set of template files and 
      their parameter correction equivalents.
      ---------------------------------------------------------------
      */
      else if(strncmp(line, "Template", strlen("Template")) == 0)
      {
         if((strstr(line, ";") == NULL) && (strstr(line, "\t") == NULL))
         {
            LogError(ERR_FILE_IO, "ParameterCorrection::CTOR(): missing separator (;) in file pair.");
         }/* end if() */

         /*--------------------------------------------------------
         Read in file pairs, taking care to preserve full path, 
         even in the presence of long and space-separated filenames.
         --------------------------------------------------------*/
         i  = ExtractColString(line, tmp1, ' '); //Template keyword
         i += ExtractFileName(&(line[i]), tmp2); //template file
         i += ExtractFileName(&(line[i]), tmp3); //parameter correction input file

         NEW_PRINT("FilePair", 1);
         pFilePair = new FilePair(tmp2, tmp3);
         MEM_CHECK(pFilePair);
         AddFilePair(pFilePair);
      }/* end if() */
   }while(strcmp(line, "EndParameterCorrection") != 0);

   fclose(pInFile);

   m_pParamGroup = pGroup;

   NEW_PRINT("ResponseVarGroup", 1);
   m_pCorrections = new ResponseVarGroup("Corrections");
   MEM_CHECK(m_pCorrections);

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
Free up memory.
******************************************************************************/
void ParameterCorrection::Destroy(void)
{
   delete m_pCorrections;
   delete m_FileList;
   delete [] m_ExecCmd;
   IncDtorCount();
}/* end Destroy() */
   
/******************************************************************************
SetExecCmd()
   Sets the syntax which is used to execute the parameter correction program
*******************************************************************************/
void ParameterCorrection::SetExecCmd(IroncladString cmd)
{
   int len;
   len = (int)strlen(cmd) + 1;
   NEW_PRINT("char", len);
   m_ExecCmd = new char[len];
   MEM_CHECK(m_ExecCmd);

   strcpy(m_ExecCmd, cmd);
} /* end SetExecCmd() */

/******************************************************************************
AddFilePair()
   Adds a file pair to the parameter correction file pair list.
******************************************************************************/
void ParameterCorrection::AddFilePair(FilePair * pFilePair)
{
   if(m_FileList == NULL) { m_FileList = pFilePair;}
   else{m_FileList->InsertPair(pFilePair);}
} /* end AddFilePair() */

/*****************************************************************************
Execute()
   Executes the external parameter correction program.
******************************************************************************/
void ParameterCorrection::Execute(void)
{  
   IroncladString dirName = GetExeDirName();
   FilePair * pCur;
   FilePipe * pPipe;

   //make substitution of parameters into parameter correction input file(s)
   pCur = m_FileList;
   while(pCur != NULL)
   {
      pPipe = pCur->GetPipe();
      m_pParamGroup->SubIntoFile(pPipe);
      pCur = pCur->GetNext();
   } /* end while() */

   //invoke system command to execute the parameter correction program
   system(m_ExecCmd);

   //extract corrected parameter values from output file(s)
   if(m_pCorrections != NULL){ m_pCorrections->ExtractVals();}
   
   //make corrections to parameter group
   double a, b;
   RespVarABC * pResp;
   ParameterABC * pParam;
   for(int i = 0; i < m_pCorrections->GetNumRespVars(); i++)
   {
      pResp = m_pCorrections->GetRespVarPtr(i);
      if(pResp != NULL)
      {
         pParam = GetParameterByName(pResp->GetName());
         if(pParam != NULL)
         {
            a = pParam->GetEstVal();
            b = pResp->GetCurrentVal();
            if(NearlyEqual(a, b) == false)
            {
	       //printf("%s %E --> %E\n", pResp->GetName(), a, b);
               pParam->SetEstVal(b);
               m_NumCorrections++;
            }/* end if() */
         }/* end if() */
      }/* end if() */
   }/* end for() */
} /* end Execute() */

/******************************************************************************
NearlyEqual()
   Test if two numbers are nearly equal to each other.
******************************************************************************/
bool ParameterCorrection::NearlyEqual(double a, double b)
{
   double absTol = 1E-6;
   double relTol = 1E-4;
   double absMax = MyMax(fabs(a), fabs(b));
   if(a == b) return true;
   double absDiff = fabs(a - b);
   if(absDiff <= absTol) return true;
   double relDiff = absDiff/absMax;
   if(relDiff <= relTol) return true;
   return false;
}/* end NearlyEqual() */

/******************************************************************************
WriteMetrics()
******************************************************************************/
void ParameterCorrection::WriteMetrics(FILE * pFile)
{  
   fprintf(pFile, "Total Parameter Corrections : %d\n", m_NumCorrections);
} /* end WriteMetrics() */

