/******************************************************************************
File     : SurrogateModel.cpp
Author   : L. Shawn Matott
Copyright: 2006, L. Shawn Matott

The SurrogateModel class encapsulates the interaction of the Ostrich 
optimization tools with externally executed groundwater modeling programs that
are simpler versions of the standard (complex) model.

Version History
04-04-06    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <string.h>
#include <stdlib.h>

#include "Model.h"
#include "SurrogateParameterGroup.h"
#include "ObservationGroup.h"
#include "ObjectiveFunction.h"
#include "FilePair.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
GetObjFuncStr()
******************************************************************************/
UnchangeableString SurrogateModel::GetObjFuncStr(void)
{
	return m_pObjFunc->GetObjFuncStr();
}/* end GetObjFuncStr() */

/******************************************************************************
 copy CTOR, copy data from an existing ModelABC and read in the rest from
 the given input file
******************************************************************************/
SurrogateModel::SurrogateModel
(
   UnmoveableString pFileName, 
   ModelABC * pComplex, 
   char * pType
)
{
   FilePair * pFilePair;
   FILE * pFile;
   char * line;
   int i, j;
   bool quoteWrap;
   char tmp1[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char tmp3[DEF_STR_SZ];
   UnmoveableString pDir = GetExeDirName();
   
   m_Counter = 0;
   m_pObsGroup = NULL;
   m_pParamGroup = NULL;
   m_pObjFunc = NULL;
   m_FileList = NULL;
   m_ExecCmd = NULL;
   m_pTypeStr = NULL;

   m_pTypeStr = new char[strlen(pType) + 1];
   strcpy(m_pTypeStr, pType);

   pFile = fopen(pFileName, "r");

   FindToken(pFile, "BeginObservations", pFileName);
   FindToken(pFile, "EndObservations", pFileName);
   rewind(pFile);
   FindToken(pFile, "BeginTiedParams", pFileName);
   FindToken(pFile, "EndTiedParams", pFileName);
   rewind(pFile);
   FindToken(pFile, "BeginFilePairs", pFileName);
   FindToken(pFile, "EndFilePairs", pFileName);
   rewind(pFile);
   FindToken(pFile, "ModelExecutable", pFileName);
   line = GetCurDataLine();

   /*--------------------------------------------------------
   Read in executable, taking care to preserve full path, 
   even in the presence of long and space-separated filenames.
   --------------------------------------------------------*/   
   i = ExtractString(line, tmp2);
   i = ValidateExtraction(i, 1, 1, "SurrogateModel");
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
   }

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
      sprintf(tmp1, "Model executable (|%s|) not found", tmp2);
      LogError(ERR_FILE_IO, tmp1);
      ExitProgram(1);
   }

   #ifdef WIN32 //windows version
      strcat(tmp1, " > ");
      strcat(tmp1, GetOstExeOut());
   #else
      strcpy(tmp2, tmp1);
      // '>&' redircts both output and error
      strcat(tmp2, " > ");
      strcat(tmp2, GetOstExeOut());
      strcat(tmp2, " 2>&1"); 
      strcpy(tmp1, tmp2);
   #endif
 
   SetCmdToExecModel(tmp1);

   /*
   --------------------------------------------------------
   Read in the 'File Pairs': a set of template files and 
   their model equivalents.
   --------------------------------------------------------
   */
   rewind(pFile);
   FindToken(pFile, "BeginFilePairs", pFileName);
   line = GetNxtDataLine(pFile, pFileName);
   while(strstr(line, "EndFilePairs") == NULL)
   {      
      if((strstr(line, ";") == NULL) && (strstr(line, "\t") == NULL))
      {
         LogError(ERR_FILE_IO, "Model::CTOR(): missing separator (;) in file pair.");
      }/* end if() */

      /*--------------------------------------------------------
      Read in file pairs, taking care to preserve full path, 
      even in the presence of long and space-separated filenames.
      --------------------------------------------------------*/
      i = ExtractFileName(line, tmp1);
      i = ExtractFileName(&(line[i]), tmp2);

      if(pDir[0] != '.')
      {
         strcpy(tmp3, pDir);
         #ifdef WIN32
            strcat(tmp3, "\\");
         #else
            strcat(tmp3, "/");
         #endif
         strcat(tmp3, tmp2);
         strcpy(tmp2, tmp3);
      }/* end if() */
      NEW_PRINT("FilePair", 1);
      pFilePair = new FilePair(tmp1, tmp2);
      MEM_CHECK(pFilePair);

      AddFilePair(pFilePair);

      line = GetNxtDataLine(pFile, pFileName);
   }/* end while() */

   fclose(pFile);

   //create a surrogate parameter group to store the tied parameters
   NEW_PRINT("SurrogateParameterGroup", 1);
   m_pParamGroup = new SurrogateParameterGroup(pFileName, 
                                               pComplex->GetParamGroupPtr());
   MEM_CHECK(m_pParamGroup);

   if(pComplex->GetObjFuncId() != OBJ_FUNC_WSSE)
   {
      LogError(ERR_IN_PARSE, "Surrogate-based calibration require WSSE objective");
      ExitProgram(1);
   }
   
   /*TO DO: need to adjust observation group to handle OST_NULL terms */
   NEW_PRINT("ObservationGroup", 1);
   m_pObsGroup = new ObservationGroup(pComplex->GetObsGroupPtr(), pFileName);
   MEM_CHECK(m_pObsGroup);
   
   //setup the objective function
   NEW_PRINT("WSSE", 1);
   m_pObjFunc = new WSSE(m_pObsGroup, false, 1.00);
   MEM_CHECK(m_pObjFunc);

   /*-----------------------------------------------------------------------
   Check template files against parameters, each parameter should appear in
   at least one template file.
   ------------------------------------------------------------------------*/
   m_pParamGroup->CheckTemplateFiles(m_FileList);

   /*-----------------------------------------------------------------------
   Delete output file, if it exists.
   ------------------------------------------------------------------------*/
   MPI_Comm_rank(MPI_COMM_WORLD, &i);
   sprintf(tmp1, "Ost%s%d.txt", m_pTypeStr, i);
   pFile = fopen(tmp1, "r");
   if(pFile != NULL)
   {
      fclose(pFile);
      if(remove(tmp1) != 0)
      {
         LogError(ERR_FILE_IO, "CTOR(): Couldn't delete file");
         ExitProgram(1);
      }
   }
   //write out banner.
   pFile = fopen(tmp1, "a+");
   if(pFile == NULL)
   {
      LogError(ERR_FILE_IO, "CTOR(): Couldn't open file");
      ExitProgram(1);
   }
   fprintf(pFile,"Run   obj. function  ");
   m_pParamGroup->Write(pFile, WRITE_BNR);
   fprintf(pFile,"\n");
   fclose(pFile);

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
Free up memory.
******************************************************************************/
void SurrogateModel::Destroy(void)
{
   delete m_pObsGroup;
   delete m_pParamGroup;
   delete m_pObjFunc;
   delete m_FileList;

   delete [] m_ExecCmd;
   delete [] m_pTypeStr;

   IncDtorCount();
}/* end Destroy() */
   
/******************************************************************************
GetObjFuncPtr() 
   Returns a pointer to the objective function.
******************************************************************************/
ObjectiveFunction * SurrogateModel::GetObjFuncPtr(void)
{
  return m_pObjFunc;
} /* end GetObjFuncPtr() */

/******************************************************************************
GetCounter()
   Returns the number of times the model has been executed
******************************************************************************/
int SurrogateModel::GetCounter(void)
{
  return m_Counter;
} /* end GetCounter() */

/******************************************************************************
SetCmdToExecModel()
   Sets the syntax which is used to execute the model
*******************************************************************************/
void SurrogateModel::SetCmdToExecModel(IroncladString cmd)
{
   int len;
   len = (int)strlen(cmd) + 1;
   NEW_PRINT("char", len);
   m_ExecCmd = new char[len];
   MEM_CHECK(m_ExecCmd);

   strcpy(m_ExecCmd, cmd);
} /* end SetCmdToExecModel() */

/******************************************************************************
AddFilePair()
   Adds a file pair to the model file pair list.
******************************************************************************/
void SurrogateModel::AddFilePair(FilePair * pFilePair)
{
   if(m_FileList == NULL) { m_FileList = pFilePair;}
   else{m_FileList->InsertPair(pFilePair);}
} /* end AddFilePair() */

/******************************************************************************
GetObsGroupPtr()
   Returns the observation group pointer
******************************************************************************/
ObservationGroup * SurrogateModel::GetObsGroupPtr(void)
{
  return m_pObsGroup;
} /* end GetObsGroupPtr() */

/*****************************************************************************
GetSurrogateParamGroupPtr()
   Returns the parameter group pointer.
******************************************************************************/
SurrogateParameterGroup * SurrogateModel::GetSurrogateParamGroupPtr(void)
{
  return m_pParamGroup;
} /* end GetSurrogateParamGroupPtr() */

/*****************************************************************************
Execute()
   Executes the surrogate model and returns the objective function value.
******************************************************************************/
double SurrogateModel::Execute(void)
{  
   IroncladString dirName = GetExeDirName();
   FilePair * pCur;
   FilePipe * pPipe;
   double val;

   //exit early if the user has requested program termination
   if(IsQuit() == true){ return NEARLY_HUGE;}

   //inc. number of times model has been executed
   m_Counter++;

   //make substitution of parameters into model input file
   pCur = m_FileList;
   while(pCur != NULL)
   {
      pPipe = pCur->GetPipe();
      m_pParamGroup->SubIntoFile(pPipe);
      pCur = pCur->GetNext();
   } /* end while() */

   //cd to model subdirectory, if needed
   if(dirName[0] != '.') { MY_CHDIR(dirName);}   

   //invoke system command to execute the model   
   system(m_ExecCmd);

   //extract computed observations from model output file(s)
   if(m_pObsGroup != NULL){ m_pObsGroup->ExtractVals();}

   //compute obj. func.
   val = m_pObjFunc->CalcObjFunc();

   //cd out of model subdirectory, if needed
   if(dirName[0] != '.') { MY_CHDIR("..");}
   
   //ouput results
   Write(val);

   m_CurObjFuncVal = val;

   return (val);
} /* end Execute() */

/******************************************************************************
Write()
   Store parameter and objective function value to model output file.
******************************************************************************/
void SurrogateModel::Write(double objFuncVal)
{
   FILE * pFile;
   char name[DEF_STR_SZ];
   int id;   

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(name, "Ost%s%d.txt", m_pTypeStr, id);

   pFile = fopen(name, "a+");
	fprintf(pFile, "%-4d  %E  ", m_Counter, objFuncVal);
   m_pParamGroup->Write(pFile, WRITE_SCI);
   fprintf(pFile, "\n");
   fclose(pFile);
} /* end Write() */

/******************************************************************************
WriteMetrics()
******************************************************************************/
void SurrogateModel::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "Total %s Evals      : %d\n", m_pTypeStr, m_Counter);
} /* end WriteMetrics() */
