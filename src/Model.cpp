/******************************************************************************
File     : Model.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The Model class encapsulates the interaction of the Ostrich optimization tools
with the externally executed modeling program. The class divides 
model components into three groups: the parameter group, the observation group 
and the objective function group. In addition to being able to execute the 
model, the Model class provides Ostrich algorithms with access to 
these groups. 

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-25-03    lsm   Modified to support long filenames in model executable.
03-05-04    lsm   added parallel model execution,
03-24-04    lsm   added support for internal models
07-08-04    lsm   added PATO hooks
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added template file check, added support for PATO.
12-02-04    lsm   Added support for Geometry parameters
01-18-05    lsm   Addded support for Generalized Constrained Optimization (GCOP)
01-18-05    lsm   Addded check for existence of model executable
10-22-05    lsm   Added support for SGI/IRIX environment, refined error handling
01-01-07    lsm   Added support for checking global sensitivity of parameters and
                  observations; items that are not sensitive over the entire range
                  of parameter values are removed from consideration. To select this
                  option, users enter the following line in the main configuration
                  section:
                     CheckSensitivities   yes
01-01-07    lsm   Added support for surrogate-model approach to calibration, ranking and
                  selection. To select this option, users enter the following line in the 
                  main configuration section:
                     SurrogateApproach   yes
01-01-07    lsm   Added model initialization and bookkepping calls.
07-13-07    lsm   Added support for the EPA SuperMUSE cluster. To select this option, 
                  users enter the following line in the main configuration section:
                     SuperMUSE yes
                  Users must also define various parameters in the 'SuperMUSE' 
                  section of the input file (See SuperMUSE.h for details):
                     BeginSuperMUSE
                       OstrichTaskerHostName <host_name>
                       TaskFile         <file_name>
                       TempFile         <file_name>
                       SuccessFile      <file_name>
                       ErrorFile        <file_name>
                       ScriptFile       <file_name>
                       MaxJobTime       <integer>
                     EndSuperMUSE
03-20-2010   lsm   Added support for warm starting Ostrich using previously
                   generated OstModel0.txt file.
******************************************************************************/
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <string>
#ifdef GCC5X
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include "Model.h"
#include "ObservationGroup.h"
#include "Observation.h"
#include "ResponseVarGroup.h"
#include "SurrogateParameterGroup.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "FilePair.h"
#include "FileList.h"
#include "AccessConverter.h"
#include "NetCDFConverter.h"
#include "PumpAndTreat.h"
#include "DecisionModule.h"
#include "SuperMUSE.h"
#include "ParameterCorrection.h"
#include "ParamInitializerABC.h"
#include "GenConstrainedOpt.h"

#include "IsoParse.h"
#include "BoxCoxModel.h"
#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"
#include "SuperMuseUtility.h"

#define JOB_SUCCEDDED (0)
#define JOB_FAILED    (1)
#define JOB_TIMED_OUT (2)

//objective function categories passed to PreserveModel output script
IroncladString ObjFuncBest          = "best";
IroncladString ObjFuncBehavioral    = "behavioral";
IroncladString ObjFuncNonBehavioral = "non-behavioral";
IroncladString ObjFuncDominated     = "dominated";
IroncladString ObjFuncNonDominated  = "non-dominated";
IroncladString ObjFuncOther         = "other";

/******************************************************************************
GetObjFuncStr()
******************************************************************************/
UnchangeableString Model::GetObjFuncStr(void)
{
	return m_pObjFunc->GetObjFuncStr(); 
}/* end GetObjFuncStr() */

/******************************************************************************
PerformParameterCorrections()
******************************************************************************/
void Model::PerformParameterCorrections(void)
{ 
	if (m_pParameterCorrection != NULL) m_pParameterCorrection->Execute(); 
}/* end PerformParameterCorrections() */

/******************************************************************************
 SaveBest()
******************************************************************************/
void Model::SaveBest(int id)
{
   char saveDir[DEF_STR_SZ];
   if(m_bSave)
   {
      //cd to model subdirectory, if needed
      if(m_DirPrefix[0] != '.') 
      {
         sprintf(saveDir, "%s%d", m_DirPrefix, id);
         MY_CHDIR(saveDir);
      }
 
      system(m_SaveCmd);

      if(m_DirPrefix[0] != '.') 
      {
         //cd out of model subdirectory, if needed
         MY_CHDIR("..");
      }
   }
}/* end SaveBest() */

/******************************************************************************
 PreserveModel()
******************************************************************************/
void Model::PreserveModel(int rank, int trial, int counter, IroncladString ofcat)
{
   char tmp[DEF_STR_SZ];

   if(m_bPreserveModelOutput == false)
   {
      return;
   }

   /* use built in preservation option
      Saves everything that was not an input file or folder
   */
   if(m_PreserveCmd == NULL)
   {   
      // Get the current path to the worker
      fs::path directory = fs::current_path();
      
      // Move up one directory to project directory
      directory = directory.parent_path();

      // Construct the archive folder
      directory /= std::string("archive");
      if (!fs::exists(directory)) {
          fs::create_directory(directory);
      };

      // Construct the worker folder
      directory /= m_DirPrefix + std::to_string(rank);
      if (!fs::exists(directory)) {
          fs::create_directory(directory);
      }

      // Construct the run folder 
      directory /= "run_" + std::to_string(counter - 1); // This adjusts by one to align with the model output files
      fs::create_directory(directory);

      // Move the files from the worker to the archive
      // TODO: Add an exclude operation to prevent copying unnecessary files
      fs::copy(fs::current_path(), directory, fs::copy_options::recursive);

      // Cleanup the worker
      m_pFileCleanupList->Cleanup(directory);

   }/* end if() */
   /* run the user-supplied preservation command */
   else
   {
      sprintf(tmp, "%s %d %d %d %s", m_PreserveCmd, rank, trial, counter, ofcat); 

      #ifdef _WIN32 //windows version
         strcat(tmp, " > OstPreserveModelOut.txt");
      #else //Linux (bash, dash, and csh)
         // '>&' redircts both output and error
         strcat(tmp, " > OstPreserveModelOut.txt 2>&1"); 
      #endif

      system(tmp);
   }
}/* end PreserveModel() */

/******************************************************************************
 default CTOR
******************************************************************************/
Model::Model(void)
{
   FilePair * pFilePair;
   char * line;
   char tmp1[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char tmp3[DEF_STR_SZ];   
   FILE * pInFile;
   int id;
   int i;
   int j;
   bool quoteWrap; //if true, then executable must be wrapped in quotes

   RegisterModelPtr(this);

   IroncladString inFileName = GetInFileName();
   UnmoveableString pDirName = GetExeDirName();

   strcpy(m_DirPrefix, ".");
   m_pParameterCorrection = NULL;
   m_pParamInitializer = NULL;
   m_pDecision = NULL;
   m_ExecCmd = NULL;
   m_SaveCmd = NULL;
   m_PreserveCmd = NULL;
   m_bSave = false;
   m_pParamGroup = NULL;
   m_pObsGroup = NULL;
   m_FileList = NULL;
   m_DbaseList = NULL;
   m_pFileCleanupList = NULL;
   m_Counter = 0;
   m_Precision = 6;
   m_pObjFunc = NULL;
   m_ObjFuncId = OBJ_FUNC_WSSE;
   m_bMultiObjProblem = false;
   m_CurMultiObjF = NULL;
   m_bCheckGlobalSens = false;
   m_bPreserveModelOutput = false;
   m_bUseSurrogates = false;
   m_CurObjFuncVal = 0.00;
   m_bWarmStart = false;
   m_bCaching = false;
   m_bDiskless = false;
   m_NumCacheHits = 0;
   m_firstCall = true;

   #ifdef _WIN32
      m_pFileCleanupList = new FileList("Ostrich.exe");
   #else
      m_pFileCleanupList = new FileList("Ostrich");
   #endif
   
   pInFile = fopen(inFileName, "r");
   if(pInFile == NULL)
   {
      FileOpenFailure("Model::CTOR", inFileName);
   }/* end if() */

   //check for critical entries, entries which have no reasonable defaults
   FindToken(pInFile, "BeginFilePairs", inFileName);
   FindToken(pInFile, "EndFilePairs", inFileName);
   rewind(pInFile);
   FindToken(pInFile, "ModelExecutable", inFileName);   

   /*
   -------------------------------------------------------------
   Read in and create the directory from which the model will be 
   run.
   -------------------------------------------------------------
   */
   rewind(pInFile);
   if((CheckToken(pInFile, "ModelSubdir", inFileName) == true) ||
      (CheckToken(pInFile, "ModelSubDir", inFileName) == true))
   {  
       line = GetCurDataLine(); 
       MyTrim(line);       
       if(strlen(line) < 12)
       { 
          LogError(ERR_IN_PARSE, "Bad ModelSubdir");
          ExitProgram(1);
       }       
       strcpy(tmp1, &line[11]);
       //strip whitespace
       MyTrim(tmp1);
       //strip quotes
       if(tmp1[0] == '"'){ tmp1[0] = ' ';}
       if(tmp1[strlen(tmp1)-1] == '"'){ tmp1[strlen(tmp1)-1] = ' ';}
       MyTrim(tmp1);
       strcpy(pDirName, tmp1);
       strcpy(m_DirPrefix, tmp1);
   }

   if(pDirName[0] != '.')
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      sprintf(tmp1, "%d", id);
      strcat(pDirName, tmp1);      
      sprintf(tmp1, "mkdir %s", pDirName);
      system(tmp1);      
   }/* end if() */

   /*
   -------------------------------------------
   Read in the model executable, modifying so 
   that output is redirected to a file. Also, 
   copy the executable file to the directory 
   of execution.
   -------------------------------------------
   */
   rewind(pInFile);
   FindToken(pInFile, "ModelExecutable", inFileName);
   line = GetCurDataLine();
   //line format = 'ModelExecutable   <var>'
   /*--------------------------------------------------------
   Read in executable, taking care to preserve full path, 
   even in the presence of long and space-separated filenames.
   --------------------------------------------------------*/   
   i = ExtractString(line, tmp2);
   i = ValidateExtraction(i, 1, 1, "Model()");
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

   /*-------------------------------------------------------------------------------------------
   Check to see if the model is internal. If it is, it will be a function call and end in a '()'
   --------------------------------------------------------------------------------------------*/
   if((strlen(tmp1) > 1) && (tmp1[strlen(tmp1)-2] == '(' ) && (tmp1[strlen(tmp1)-1] == ')' ))
   {
      m_InternalModel = true;
   }
   else
   {
      m_InternalModel = false;
   }

   /*------------------------------------------------------------------------------------------
   Extract executable file name from rest of path and save it. This is so that Ostrich can 
   cleanup after itself when copying files around.
   ------------------------------------------------------------------------------------------*/
   if(m_InternalModel == false)
   {
      j = 0;
      for(i = (int)(strlen(tmp1))-1; i > 0; i--)
      {
         if((tmp1[i] != '\\') && (tmp1[i] != '/'))
         {
            tmp2[j] = tmp1[i];
            j++;
         }
         else
         {
            if(quoteWrap == true)
            {
               tmp2[j] = '"';
               j++;
            }
            tmp2[j] = (char)NULL;
            MyStrRev(tmp2);
            m_pFileCleanupList->Insert(tmp2);
            break;
         }
      }
   }/* end if() */

   if((pDirName[0] != '.') && (m_InternalModel == false))
   {      
      #ifdef _WIN32
         sprintf(tmp2, "copy %s %s", tmp1, pDirName);
      #else
         sprintf(tmp2, "cp %s %s", tmp1, pDirName);
      #endif
      system(tmp2);
   } /* end if() */

   if(m_InternalModel == false)
   {
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
         LogError(ERR_FILE_IO, tmp1);
         ExitProgram(1);
      }

      #ifdef _WIN32 //windows version
         strcat(tmp1, " > ");
         strcat(tmp1, GetOstExeOut());
      #else //Linux (bash, dash, csh)
         strcpy(tmp2, tmp1);
         // '>&' redircts both output and error
         strcat(tmp2, " > ");
         strcat(tmp2, GetOstExeOut());
         strcat(tmp2, " 2>&1"); 
         strcpy(tmp1, tmp2);
      #endif
   }
   SetCmdToExecModel(tmp1);

   /*
   --------------------------------------------------------
   Read in the 'File Pairs': a set of template files and 
   their model equivalents.
   --------------------------------------------------------
   */
   rewind(pInFile);
   FindToken(pInFile, "BeginFilePairs", inFileName);
   line = GetNxtDataLine(pInFile, inFileName);
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
      //first file in pair.....
      i = ExtractFileName(line, tmp1);
      i = ExtractFileName(&(line[i]), tmp2);

      if(pDirName[0] != '.')
      {
         strcpy(tmp3, pDirName);
         #ifdef _WIN32
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

      line = GetNxtDataLine(pInFile, inFileName);
   }/* end while() */

   /*
   --------------------------------------------------------------------
   Read in any extra model files, these will need to be copied to the 
   model subdirectory.
   --------------------------------------------------------------------
   */
   rewind(pInFile);
   if(CheckToken(pInFile, "BeginExtraFiles", inFileName) == true)
   {
      //make sure end token exists
      FindToken(pInFile, "EndExtraFiles", inFileName);
      rewind(pInFile);
      FindToken(pInFile, "BeginExtraFiles", inFileName);

      line = GetNxtDataLine(pInFile, inFileName);
      while(strstr(line, "EndExtraFiles") == NULL)
      {
         //extra file
         ExtractFileName(line, tmp1);

         // add to cleanup list
         m_pFileCleanupList->Insert(tmp1);

         // Move the file
         if(pDirName[0] != '.') {
             // Construct the worker path 
             fs::path workerPath = fs::current_path() /= std::string(pDirName);
             workerPath /= tmp1;

             // Construct the path to the source file
             fs::path sourcePath = fs::current_path() /= tmp1;

            // Copy the file from the source to the worker path
             fs::copy_file(sourcePath, workerPath);
               
         }/* end if() */

         line = GetNxtDataLine(pInFile, inFileName);
      }/* end while() */
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in any extra model folders, these will need to be copied (along
   with their contects) to the model subdirectory.
   --------------------------------------------------------------------
   */
   rewind(pInFile);
   if(CheckToken(pInFile, "BeginExtraDirs", inFileName) == true)
   {
      //make sure end token exists
      FindToken(pInFile, "EndExtraDirs", inFileName);
      rewind(pInFile);
      FindToken(pInFile, "BeginExtraDirs", inFileName);

      line = GetNxtDataLine(pInFile, inFileName);
      while(strcmp(line, "EndExtraDirs") != 0)
      {
         //extra dir
         ExtractFileName(line, tmp1);

         // Move all entries in the directory
         if(pDirName[0] != '.') {
             // Construct the worker path 
             fs::path workerPath = fs::current_path() /= std::string(pDirName);
             workerPath /= tmp1;

             if (~fs::exists(workerPath)) {
                 fs::create_directories(workerPath);
             }

             // Construct the path to the source source
             fs::path sourcePath = fs::current_path() /= tmp1;

             // Copy the files
             fs::copy(sourcePath, workerPath, fs::copy_options::recursive);

             // Get a list of files in the source directory
             std::vector<std::string> sourceFiles;
             for (const auto& entry : fs::recursive_directory_iterator(sourcePath)) {
                 if (!fs::is_directory(entry.path())) {
                     fs::path partialPath = entry.path().lexically_relative(fs::current_path());
                     sourceFiles.push_back(partialPath.string());
                 }
             }

             // Add files within the directory to the cleanup list
             for (int entry=0; entry < sourceFiles.size(); entry++) {
                 m_pFileCleanupList->Insert(&sourceFiles[entry][0]);             
             }

         }/* end if() */

         line = GetNxtDataLine(pInFile, inFileName);
      }/* end while() */
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in diskless model
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "DisklessModel", inFileName) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) {m_bDiskless = true;}
   }/* end if() */

   /*
   -------------------------------------------
   Read in the name of the batch file to run
   whenever a new 'best' solution is found.
   -------------------------------------------
   */
   rewind(pInFile);
   if(CheckToken(pInFile, "PreserveBestModel", inFileName) == true)
   {
      line = GetCurDataLine();
      //line format = 'PreserveBestModel   <var>'
      /*--------------------------------------------------------
      Read in executable, taking care to preserve full path, 
      even in the presence of long and space-separated filenames.
      --------------------------------------------------------*/   
      i = ExtractString(line, tmp2);
      i = ValidateExtraction(i, 1, 1, "Model()");
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

      if(pDirName[0] != '.')
      {
         #ifdef _WIN32
            sprintf(tmp2, "copy %s %s", tmp1, pDirName);
         #else
            sprintf(tmp2, "cp %s %s", tmp1, pDirName);
         #endif
         system(tmp2);
      } /* end if() */

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
         sprintf(tmp1, "File for saving best solution (|%s|) not found", tmp2);
         LogError(ERR_FILE_IO, tmp1);
         ExitProgram(1);
      }

      #ifdef _WIN32 //windows version
         strcat(tmp1, " > OstSaveOut.txt");
      #else //Linux (bash, dash, and csh)
         strcpy(tmp2, tmp1);
         // '>&' redircts both output and error
         strcat(tmp2, " > OstSaveOut.txt 2>&1"); 
         strcpy(tmp1, tmp2);
      #endif
      m_bSave = true;
      m_SaveCmd = new char[(int)strlen(tmp1) + 1];
      strcpy(m_SaveCmd, tmp1);
   }/* end if(PreserveBestModel) */

   /*
   --------------------------------------------------------------------
   Check for alternate objective function (default is WSSE)
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "ObjectiveFunction", inFileName) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strstr(tmp2, "user") != NULL) {m_ObjFuncId = OBJ_FUNC_USER;}
      else if(strstr(tmp2, "sawe") != NULL) {m_ObjFuncId = OBJ_FUNC_SAWE;}
      else if(strstr(tmp2, "wsse") != NULL) {m_ObjFuncId = OBJ_FUNC_WSSE;}
      else if(strstr(tmp2, "pato") != NULL) {m_ObjFuncId = OBJ_FUNC_PATO;}
      else if(strstr(tmp2, "gcop") != NULL) {m_ObjFuncId = OBJ_FUNC_GCOP;}
   } /* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to check sensitivities
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "CheckSensitivities", inFileName) == true)
   {   
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) {m_bCheckGlobalSens = true;}
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to use surrogate models
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "SurrogateApproach", inFileName) == true)
   {   
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) {m_bUseSurrogates = true;}
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to attempt file cleanup - seems to be rather buggy at
   the moment so lets give users the ability to opt out
   --------------------------------------------------------------------
   */
   m_bTryFileCleanup = false;
   rewind(pInFile);
   if(CheckToken(pInFile, "AttemptFileCleanup", inFileName) == true)
   {   
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) {m_bTryFileCleanup = true;}
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to use SuperMUSE
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   bool bSMUSE = false;
   if(CheckToken(pInFile, "SuperMUSE", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         EnableSuperMUSE();
         InitSuperMUSE(pInFile, (ModelABC *)this);
         bSMUSE = true;
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to preserve model output files. This only applies if
   a model subdirectory is used.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "PreserveModelOutput", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      MyTrim(tmp2);
      if(strcmp(tmp2, "yes") == 0) 
      {
         m_bPreserveModelOutput = true;
      }
      else if(strcmp(tmp2, "no") == 0) 
      {
         m_bPreserveModelOutput = false;
      }
      else //a valid PreserveModelOutput script?
      {
         //line format = 'PreserveModelOutput  <var>'
         /*--------------------------------------------------------
         Read in executable, taking care to preserve full path, 
         even in the presence of long and space-separated filenames.
         --------------------------------------------------------*/   
         i = ExtractString(line, tmp2);
         i = ValidateExtraction(i, 1, 1, "Model()");
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

         //stage to workdir if needed
         if(pDirName[0] != '.')
         {
            #ifdef _WIN32
               sprintf(tmp2, "copy %s %s", tmp1, pDirName);
            #else
               sprintf(tmp2, "cp %s %s", tmp1, pDirName);
            #endif
            system(tmp2);
         } /* end if() */

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
            sprintf(tmp1, "File for preserving model output (|%s|) not found", tmp2);
            LogError(ERR_FILE_IO, tmp1);
            ExitProgram(1);
         }

         m_bPreserveModelOutput = true;
         m_PreserveCmd = new char[(int)strlen(tmp1) + 1];
         strcpy(m_PreserveCmd, tmp1);
      } /* end else() */
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in warm start flag. This only applies if we want to restart
   Ostrich from a previously aborted or interrupted run in which
   the OstModel0.txt file has been preserved.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "OstrichWarmStart", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         printf("Warm Start has been activated\n");
         printf("Ostrich will resume a pervious search.\n");
         m_bWarmStart = true;
         RestoreRandomSeed();
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in caching flag. This only applies if we want to search the
   history of model evaluations (stored in OstModel0.txt) prior to
   running the model. If the candidate solution has already been
   evaluated, it will be found in the cache and so we won't have to
   run the model.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "OstrichCaching", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         m_bCaching = true;
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in number of digits of precision in I/O.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "NumDigitsOfPrecision", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %d", tmp1, &m_Precision);
      if((m_Precision < 1) || (m_Precision > 32)) 
      {
         LogError(ERR_FILE_IO, "Invalid precision setting - defaulting to 6 digits.");
         m_Precision = 6;
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in telescoping strategy, if any.
   --------------------------------------------------------------------
   */
   char tscp[DEF_STR_SZ];
   rewind(pInFile);
   if(CheckToken(pInFile, "TelescopingStrategy", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tscp);

      MyStrLwr(tscp);
      if(strcmp(tscp, "convex-power")    == 0) m_Telescope = TSCOPE_PVEX;
      if(strcmp(tscp, "convex")          == 0) m_Telescope = TSCOPE_CVEX;
      if(strcmp(tscp, "linear")          == 0) m_Telescope = TSCOPE_LINR;
      if(strcmp(tscp, "concave")         == 0) m_Telescope = TSCOPE_CAVE;
      if(strcmp(tscp, "delayed-concave") == 0) m_Telescope = TSCOPE_DCVE;
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in BoxCox transformation, if any.
   --------------------------------------------------------------------
   */
   char boxcox[DEF_STR_SZ];
   bool bBoxCox = false;
   double boxCoxValue=1.00;
   rewind(pInFile);
   if(CheckToken(pInFile, "BoxCoxTransformation", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, boxcox);
      MyStrLwr(boxcox);
      if(strcmp(boxcox, "extract") == 0)
      {
         boxCoxValue = ExtractBoxCoxValue();
      }
      else
      {
         boxCoxValue = atof(boxcox);
      }
      bBoxCox = true;
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in Parameter Initializer type, if any.
   --------------------------------------------------------------------
   */
   char paramInitializer[DEF_STR_SZ];
   strcpy(paramInitializer, "");
   rewind(pInFile);
   if(CheckToken(pInFile, "ParamInitializer", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, paramInitializer);
      MyStrLwr(paramInitializer);
   }/* end if() */

   fclose(pInFile);

   if(bSMUSE) CleanSuperMUSE();

   /*
   --------------------------------------------------------
   Read in database conversion information.
   --------------------------------------------------------
   */
	DatabaseABC * access_DbaseList = new AccessConverter();
   DatabaseABC * netcdf_DbaseList = new NetCDFConverter();
	MEM_CHECK(access_DbaseList);
	MEM_CHECK(netcdf_DbaseList);
	if(access_DbaseList->ReadFromFile())
   {
		NEW_PRINT("AccessConverter", 1);
		m_DbaseList = access_DbaseList;
      delete netcdf_DbaseList;
		netcdf_DbaseList = NULL;
   }
	else if(netcdf_DbaseList->ReadFromFile())
	{
		NEW_PRINT("NetCDFConverter", 1);
		m_DbaseList = netcdf_DbaseList;
		delete access_DbaseList;
		access_DbaseList = NULL;
	}

   NEW_PRINT("ParameterGroup", 1);
   m_pParamGroup = new ParameterGroup();
   MEM_CHECK(m_pParamGroup);

   //check bounds on parameter group
   m_pParamGroup->CheckBounds();

   if((m_ObjFuncId == OBJ_FUNC_WSSE) || (m_ObjFuncId == OBJ_FUNC_SAWE))
   {
      NEW_PRINT("ObservationGroup", 1);
      m_pObsGroup = new ObservationGroup();
      MEM_CHECK(m_pObsGroup);
   }/* end if() */   

   //setup the objective function
   if(m_ObjFuncId == OBJ_FUNC_WSSE)
   { 
      NEW_PRINT("WSSE", 1);
      m_pObjFunc = new WSSE(m_pObsGroup, bBoxCox, boxCoxValue);
   }
   else if(m_ObjFuncId == OBJ_FUNC_SAWE)
   { 
      NEW_PRINT("SAWE", 1);
      m_pObjFunc = new SAWE(m_pObsGroup);
   }
   else if(m_ObjFuncId == OBJ_FUNC_PATO)
   { 
      NEW_PRINT("PATO", 1);
      m_pObjFunc = new PATO(m_pParamGroup);
   }
   else if(m_ObjFuncId == OBJ_FUNC_GCOP)
   { 
      NEW_PRINT("GCOP", 1);
      m_pObjFunc = new GCOP(m_pParamGroup);
   }
   else /* OBJ_FUNC_USER) */ 
   { 
      NEW_PRINT("USER", 1);
      m_pObjFunc = new UserObjFunc(GetOstExeOut());
   }
   MEM_CHECK(m_pObjFunc);

   if(m_pObjFunc->CalcMultiObjFunc(NULL, -1) > 1)
   {
      m_bMultiObjProblem = true;
   }

   /*-----------------------------------------------------------------------
   Read in special parameters.
   ------------------------------------------------------------------------*/
   m_pParamGroup->InitSpecialParams(inFileName);

   /*-----------------------------------------------------------------------
   Check template files against parameters, each parameter should appear in
   at least one template file or at least one database entry.
   ------------------------------------------------------------------------*/
   m_pParamGroup->CheckTemplateFiles(m_FileList);
   //m_pParamGroup->CheckDbaseFiles(m_DbaseList);

   /*-----------------------------------------------------------------------
   Check parameters for uniqueness, each parameter should be unique and 
   should not be a substring of another parameter.
   ------------------------------------------------------------------------*/
   m_pParamGroup->CheckMnemonics();

   /*-----------------------------------------------------------------------
   Initialize surrogate models, if the surrogate-based approach is enabled.
   ------------------------------------------------------------------------*/
   if(m_bUseSurrogates == true)
   {
      m_pDecision = new DecisionModule((ModelABC *)this);
   }

   ExcludeConstantParameters();

   CheckGlobalSensitivity();
   
   // check for parameter correction strategy
   pInFile = fopen(inFileName, "r");
   if(pInFile == NULL)
   {
      FileOpenFailure("Model::CTOR", inFileName);
   }/* end if() */
   if(CheckToken(pInFile, "BeginParameterCorrection", inFileName) == true)
   {
      m_pParameterCorrection = new ParameterCorrection(m_pParamGroup);
   }/* end if() */
   fclose(pInFile);

   // check for parameter initialization strategy
   pInFile = fopen(inFileName, "r");
   if(pInFile == NULL)
   {
      FileOpenFailure("Model::CTOR", inFileName);
   }/* end if() */
   if(strcmp(paramInitializer, "hamed") == 0)
   {
      m_pParamInitializer = new HamedParamInitializer(m_pParamGroup, pInFile);
   }
   else if(strcmp(paramInitializer, "kmeans") == 0)
   {
      m_pParamInitializer = new KMeansParamInitializer(m_pParamGroup, pInFile);
   }
   fclose(pInFile);

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
ExtractBoxCoxValue()

Retrieve optimal BoxCox transformation from previous OstOutput0.txt file.
******************************************************************************/
double ExtractBoxCoxValue(void)
{
   FILE * pOut;
   char * line;
   int max_line_size;
   const char * pTok = "Estimated Optimal Box-Cox Transformation";
   double b;

   max_line_size = GetMaxLineSizeInFile((char *)"OstOutput0.txt");
   line = new char[max_line_size + 1];
   pOut = fopen("OstOutput0.txt", "r");
   if(pOut == NULL)
   {
      LogError(ERR_FILE_IO, "Unable to extract Box-Cox transformation value. Defaulting to 1.00.");
      delete [] line;
      return 1.00;
   }
   while(!feof(pOut))
   {
      fgets(line, max_line_size, pOut);
      if(strncmp(line, pTok, strlen(pTok)) == 0)
      {
         fgets(line, max_line_size, pOut);
         sscanf(line, "Lambda : %lf\n", &b);
         fclose(pOut);
         delete [] line;
         return b;
      }
   }/* end while() */
   fclose(pOut);
   LogError(ERR_FILE_IO, "Unable to extract Box-Cox transformation value. Defaulting to 1.00.");
   delete [] line;
   return 1.00;   
}/* end ExtractBoxCoxValue() */

/******************************************************************************
Fetch the parameter initializer.
******************************************************************************/
ParamInitializerABC * Model::GetParamInitializerPtr(void)
{
   return m_pParamInitializer;
}/* end GetParamInitializerPtr() */

/******************************************************************************
Free up memory.
******************************************************************************/
void Model::Destroy(void)
{
   delete m_pObsGroup;
   delete m_pParamGroup;
   delete m_pParameterCorrection;
   delete m_pParamInitializer;
   delete m_pObjFunc;
   delete m_FileList;
   delete m_DbaseList;
   delete [] m_ExecCmd;
   delete [] m_SaveCmd;
   delete [] m_PreserveCmd;
   delete [] m_CurMultiObjF;
   m_bSave = false;
   delete m_pDecision;

   if((m_pFileCleanupList != NULL) && (m_bTryFileCleanup == true))
   {
      IroncladString dirName = GetExeDirName(); 
      if(dirName[0] != '.')
      {
         m_pFileCleanupList->Cleanup(dirName);         
      }
      delete m_pFileCleanupList;
   }

   //cleanup diskless data if needed
   DisklessIsotherm(NULL, NULL);
   DisklessMcCammon(NULL, NULL);

   IncDtorCount();
}/* end Destroy() */
   
/******************************************************************************
GetObjFuncPtr() 
   Returns a pointer to the objective function.
******************************************************************************/
ObjectiveFunction * Model::GetObjFuncPtr(void)
{
  return m_pObjFunc;
} /* end GetObjFuncPtr() */

/******************************************************************************
GetCounter()
   Returns the number of times the model has been executed
******************************************************************************/
int Model::GetCounter(void)
{
  return m_Counter;
} /* end GetCounter() */

/******************************************************************************
SetCounter()
******************************************************************************/
void Model::SetCounter(int count)
{
   m_Counter = count;
   if(m_Counter == 0)
   {
      m_firstCall = true;
   }
   else
   {
      m_firstCall = false;
   }
} /* end SetCounter() */

/******************************************************************************
SetCmdToExecModel()
   Sets the syntax which is used to execute the model
*******************************************************************************/
void Model::SetCmdToExecModel(IroncladString cmd)
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
void Model::AddFilePair(FilePair * pFilePair)
{
   if(m_FileList == NULL) { m_FileList = pFilePair;}
   else{m_FileList->InsertPair(pFilePair);}
} /* end AddFilePair() */

/******************************************************************************
AddDatabase()
   Adds a database conversion to the list.
******************************************************************************/
void Model::AddDatabase(DatabaseABC * pDbase)
{
   if(m_DbaseList == NULL) { m_DbaseList = pDbase;}
   else{m_DbaseList->InsertDbase(pDbase);}
} /* end AddDatabase() */

/******************************************************************************
GetObsGroupPtr()
   Returns the observation group pointer
******************************************************************************/
ObservationGroup * Model::GetObsGroupPtr(void)
{
  return m_pObsGroup;
} /* end GetObsGroupPtr() */

/*****************************************************************************
GetParamGroupPtr()
   Returns the parameter group pointer.
******************************************************************************/
ParameterGroup * Model::GetParamGroupPtr(void)
{
  return m_pParamGroup;
} /* end GetParamGroupPtr() */

/*****************************************************************************
Execute()
   Executes the model (or surrogate) and returns the objective function value.
******************************************************************************/
double Model::Execute(void)
{  
   double F;

   if((m_bDiskless == true) && (m_InternalModel == true) && (m_Counter > 0))
   {
      F = DisklessExecute();
   }
   else if(m_bUseSurrogates == false)
   {
      F = StdExecute(0.00);
   } 
   else
   {
      F = m_pDecision->Execute();
   }

   //if desired update log of residuals
   WriteIterationResiduals();

   return F;
} /* end Execute() */

/*****************************************************************************
Execute()
   Executes the model returns a vector of objective function values.
******************************************************************************/
void Model::Execute(double * pF, int nObj)
{  
   IroncladString dirName = GetExeDirName();
   FilePair * pCur;
   FilePipe * pPipe;
   int rank;

   //initialize costs
   for(int i = 0; i < nObj; i++)
   {
      pF[i] = NEARLY_HUGE;
   }

   //exit early if the user has requested program termination
   if(IsQuit() == true)
   {
      return;
   }

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

   //make substitution of parameters into model input databases
   if(m_DbaseList != NULL)
   {
      m_pParamGroup->SubIntoDbase(m_DbaseList);
   } /* end while() */

   //invoke system command to execute the model   
   system(m_ExecCmd);

   //extract computed reponses from model output database(s)
   if(m_DbaseList != NULL)
   {
      DatabaseABC * pCur;

      //first clean up ASCII files
      for(pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
      {
         pCur->DeleteASCIIFile();
      }

      //convert the responses to ASCII file
      for(pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
      {
         pCur->ReadResponse();
      }
   }

   //extract computed observations from model output file(s)
   if(m_pObsGroup != NULL){ m_pObsGroup->ExtractVals();}

   //compute obj. func.
   m_pObjFunc->CalcMultiObjFunc(pF, nObj);

   //categorize the obj. func.
   if (dirName[0] != '.') { MY_CHDIR(".."); }
   IroncladString pCatStr = GetObjFuncCategory(pF, nObj);
   if (dirName[0] != '.') { MY_CHDIR(dirName); }

   //preserve model output, if desired
   PreserveModel(rank, GetTrialNumber(), m_Counter, pCatStr);

   //cd out of model subdirectory, if needed
   if(dirName[0] != '.') { MY_CHDIR("..");}
   
   //ouput results (first objective, others should be tagged as augmented responses)
   double val = pF[0];
   m_CurObjFuncVal = val;

   //store copy of latest result, needed for printing purposes
   if(m_CurMultiObjF == NULL) m_CurMultiObjF = new double[nObj];
   for(int i = 0; i < nObj; i++)
   {
      m_CurMultiObjF[i] = pF[i];
   }

   Write(val);
} /* end Execute() */

/*****************************************************************************
DisklessExecute()
   Executes an inernal model without using I/O.
******************************************************************************/
double Model::DisklessExecute(void)
{   
   //inc. number of times model has been executed
   m_Counter++;

   if(strcmp(m_ExecCmd, "Isotherm()") == 0){ DisklessIsotherm(m_pParamGroup, m_pObsGroup);}
   //else if(strcmp(m_ExecCmd, "Orear()") == 0){ DisklessOrear();}
   else if(strcmp(m_ExecCmd, "McCammon()") == 0){ DisklessMcCammon(m_pParamGroup, m_pObsGroup);}
   else if(strcmp(m_ExecCmd, "Kinniburgh()") == 0){ DisklessKinniburgh(m_pParamGroup, m_pObsGroup);}
   //else if(strcmp(m_ExecCmd, "AdvancedKinniburgh()") == 0){ DisklessAdvancedKinniburgh();}
   //else if(strcmp(m_ExecCmd, "BoxCox()") == 0){ DisklessBoxCoxModel();}
   else{ return 0.00;}

   //compute obj. func.
   double val = m_pObjFunc->CalcObjFunc();
 
   m_CurObjFuncVal = val;

   return (val);
}/* end DisklessExecute() */

/*****************************************************************************
Execute()
   Executes the model (or surrogate) and returns the objective function value.
   Incorporates a penalty function for violation of parameter bounds.

   NOT compatible with surrogate-based appraoch.
******************************************************************************/
double Model::Execute(double viol)
{
	return(StdExecute(viol));
} /* end Execute() */


/*****************************************************************************
GetObjFuncCategory()
   Determine the obj. func. category argument to pass to the PreserveModel
   script.
******************************************************************************/
IroncladString Model::GetObjFuncCategory(double * pF, int nObj)
{
   switch(GetProgramType())
   {
      /* single objective optimizers */
      case(GA_PROGRAM):
      case(BGA_PROGRAM):
      case(SA_PROGRAM):
      case(CSA_PROGRAM):
      case(VSA_PROGRAM):
      case(PSO_PROGRAM):
      case(PSO_LEV_PROGRAM):
      case(LEV_PROGRAM):
      case(POWL_PROGRAM):
      case(BIS_PROGRAM):
      case(STEEP_PROGRAM):
      case(FLRV_PROGRAM):
      case(DDS_PROGRAM):
      case(GMLMS_PROGRAM): 
      case(SCEUA_PROGRAM):
      case(DDDS_PROGRAM):
      case(SMP_PROGRAM):
      case(PDDS_PROGRAM):
      case(APPSO_PROGRAM):
      case(BEERS_PROGRAM):
      {
         if((m_Counter <= 1) || (pF[0] < GetBestObjFunc(m_pParamGroup->GetNumParams())))
         {
            return(ObjFuncBest);
         }
         return ObjFuncOther;
      }/* end single obj. optimizers */

      /* rejection samplers */
      case(RJSMP_PROGRAM):
      case(METRO_PROGRAM):
      {
         return ObjFuncOther;
      }/* end rejection samplers */

      /* behaviorial samplers */
      case(GLUE_PROGRAM):
      case(DDSAU_PROGRAM):
      {
         if(pF[0] < GetObjFuncThreshold())
         {
            return(ObjFuncBehavioral);
         }
         return ObjFuncNonBehavioral;
      }/* end behavioral samplers */

      /* multi-objective optimizers */
      case(SMOOTH_PROGRAM):
      case(MOPSOCD_PROGRAM):
      case(NSGAII_PROGRAM):
      case(PAES_PROGRAM):
      case(PADDS_PROGRAM):
      case(PARA_PADDS_PROGRAM):
     {
         if((m_Counter <= 1) || IsNonDominated(pF, nObj))
         {
            return ObjFuncNonDominated;
         }
         return ObjFuncDominated;
      }/* end multi-objective optimizers*/

      /* miscellaneous*/
      case(SET_INFILE):
      case(STATS_PROGRAM):
      case(UTIL_PROGRAM):
      case(GRID_PROGRAM):
      case(EVAL_PROGRAM):
      case(JACOBIAN_PROGRAM):
      case(HESSIAN_PROGRAM):
      case(GRADIENT_PROGRAM):
      case(QUIT_PROGRAM):
      default:
      {
         return ObjFuncOther;
      }/* end miscellaneous */
   }/* end case() */
   return ObjFuncOther;
}/* end GetObjFuncCategory() */

/*****************************************************************************
Bookkeep()
   Performs bookkeeping operations related to parallel executing.
   if bFinal == true, iteration is complete so collect metrics
   if bFinal == false, in the middle of iteration, share information between
    processors
******************************************************************************/
void Model::Bookkeep(bool bFinal)
{  
   int id, nprocs, temp, i;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (nprocs == 1) return;

   if(m_bUseSurrogates == true)
   {
      m_pDecision->Bookkeep(bFinal);
   }   

   /* ---------------------------------------------------------------
   Collect total evals 
   --------------------------------------------------------------- */
   if(bFinal == true)
   {
      for(i = 1; i < nprocs; i++)
      {
         temp = m_Counter;
         MPI_Bcast(&temp, 1, MPI_INTEGER, i, MPI_COMM_WORLD);
         if(id == 0) m_Counter += temp;
      }
   }
} /* end Bookkeep() */

/*****************************************************************************
StdExecute()
   Executes the standard (complex) model and returns the objective function 
   value.
******************************************************************************/
double Model::StdExecute(double viol)
{   
   IroncladString dirName = GetExeDirName();
   FilePair * pCur;
   FilePipe * pPipe;
   double val;
   int rank;
   bool isGoodTopo;

   //exit early if the user has requested program termination
   if(IsQuit() == true){ return NEARLY_HUGE;}

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   //inc. number of times model has been executed
   m_Counter++;

   /*
   adjust geometries to conform to topology rules, 
   if the topology cannot be 'fixed' then false will 
   be returned....
   */
   isGoodTopo = m_pParamGroup->FixGeometry();
   if(isGoodTopo == false)
   {
      LogError(ERR_MODL_EXE, "Could not correct model topology");
   }

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

   //make substitution of parameters into model input databases
   if(m_DbaseList != NULL)
   {
      m_pParamGroup->SubIntoDbase(m_DbaseList);
   } /* end while() */

   /* -----------------------------------------------------------
   If caching is enabled, attempt to read model evaluation
   from OstModel0.txt file.
   ----------------------------------------------------------- */
   if(m_bCaching == true)
   {      
      bool bCached = CheckCache(&val);
      if(bCached == true) //found previous model result
      {
         m_NumCacheHits++;
         Write(val);
         m_CurObjFuncVal = val;
         return (val);
      }
   }/* end if() */

   if(m_InternalModel == true)
   {
      if(strcmp(m_ExecCmd, "Isotherm()") == 0){ Isotherm(m_bDiskless);}
      else if(strcmp(m_ExecCmd, "Orear()") == 0){ Orear();}
      else if(strcmp(m_ExecCmd, "McCammon()") == 0){ McCammon(m_bDiskless);}
      else if(strcmp(m_ExecCmd, "Kinniburgh()") == 0){ Kinniburgh(m_bDiskless);}
      else if(strcmp(m_ExecCmd, "AdvancedKinniburgh()") == 0){ AdvancedKinniburgh();}
      else if(strcmp(m_ExecCmd, "BoxCox()") == 0){ BoxCoxModel();}
      else{LogError(ERR_BAD_ARGS, "Unknown internal model"); ExitProgram(1);}
   }
   else
   {
      //invoke system command to execute the model   
      system(m_ExecCmd);
   }/* end else (external model) */

   //extract computed reponses from model output database(s)
   if(m_DbaseList != NULL)
   {
      DatabaseABC * pCur;

      //first clean up ASCII files
      for(pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
      {
         pCur->DeleteASCIIFile();
      }

      //convert the responses to ASCII file
      for(pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
      {
         pCur->ReadResponse();
      }
   }

   //extract computed observations from model output file(s)
   if(m_pObsGroup != NULL){ m_pObsGroup->ExtractVals();}

   //compute obj. func.
   val = m_pObjFunc->CalcObjFunc();

   //add in penalty for violation of parameter bounds
   val += viol*MyMax(1.00,val);

   //categorize the obj. func.
   IroncladString pCatStr = GetObjFuncCategory(&val, 1);

   //preserve model output, if desired
   PreserveModel(rank, GetTrialNumber(), m_Counter, pCatStr);

   //cd out of model subdirectory, if needed
   if(dirName[0] != '.') { MY_CHDIR("..");}
   
   //ouput results
   Write(val);

   m_CurObjFuncVal = val;

   return (val);
} /* end StdExecute() */

/******************************************************************************
CheckCache()
   Attempt to read previous result stored in OstModel0.txt file. If successful
   the corresponding obj. function value will be stored in val argument. Other-
   wise the function returns false and val is unchanged.
******************************************************************************/
bool Model::CheckCache(double * val)
{
   int i, j, max_line_size;
   char * line;
   char valStr[DEF_STR_SZ];
   char * pTmp;
   double objFunc, paramVal;
   bool bFound;
   FILE * pIn;

   ParameterABC * pParam;
   ParameterGroup * pGroup = GetParamGroupPtr();
   if(pGroup == NULL) return false;
   int np = pGroup->GetNumParams();

   max_line_size = GetMaxLineSizeInFile((char *)"OstModel0.txt");
   line = new char[max_line_size+1];
   line[0] = NULLSTR;
   pIn = fopen("OstModel0.txt", "r");
   if(pIn == NULL)
   {
      delete [] line;
      return false;
   }
   while(!feof(pIn))
   {
      fgets(line, max_line_size, pIn);
      pTmp = line;
      MyTrim(pTmp);
      j = ExtractColString(pTmp, valStr, ' ');
      pTmp += j;
      MyTrim(pTmp);
      j = ExtractColString(pTmp, valStr, ' ');
      pTmp += j;
      MyTrim(pTmp);
      objFunc = atof(valStr);

      bFound = true;
      for(i = 0; i < np; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         j = ExtractColString(pTmp, valStr, ' ');
         paramVal = pParam->ConvertInVal(atof(valStr));
         if(fabs(paramVal - pParam->GetEstVal()) > 1E-10)
         {
            bFound = false;
            break;
         }/* end if() */
         pTmp += j;
         MyTrim(pTmp);
      }/* end for() */
      if(bFound == true)
      {
         fclose(pIn);
         *val = objFunc;
         delete [] line;
         return true;
      }
   }/* end while() */
   fclose(pIn);
   delete [] line;
   return false;
}/* end CheckCache() */

/******************************************************************************
GatherTask()
   Read output file of a SuperMUSE task (stored in pDir directory) and compute 
   the associated objective function.
******************************************************************************/
double Model::GatherTask(char * pDir)
{
  double val;

  //inc. number of times model has been executed
  m_Counter++;

  //cd to task subdirectory
  MY_CHDIR(pDir);

  //extract computed observations from model output file(s)
  if(m_pObsGroup != NULL){ m_pObsGroup->ExtractVals();}

  //compute obj. func.
  val = m_pObjFunc->CalcObjFunc();

  //cd out of task subdirectory
  MY_CHDIR("..");

  //ouput results
  Write(val);

  m_CurObjFuncVal = val;

  return (val);
}/* end GatherTask() */

/******************************************************************************
Write()
   Store parameter and objective function value to model output file.
******************************************************************************/
void Model::Write(double objFuncVal)
{
   ResponseVarGroup * pRespVarGroup;
   FILE * pFile;
   char name[DEF_STR_SZ];
   int id;   

   pRespVarGroup = NULL;
   if(m_pObjFunc != NULL)
   {
      pRespVarGroup = (ResponseVarGroup *)(m_pObjFunc->GetResponseVarGroup());
   }

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(name, "OstModel%d.txt", id);

   if(m_firstCall == true) 
   {
      pFile = fopen(name, "r");
      if(pFile != NULL)
      {
         fclose(pFile);
         //don't remove if using WarmStart
         if(m_bWarmStart == false)
         {
            if(remove(name) != 0)
            {
               LogError(ERR_FILE_IO, "Write(): Couldn't delete OstModel.txt file");
               ExitProgram(1);
            }
         }
      }
      m_firstCall = false;
      //write out banner.
      pFile = fopen(name, "a+");
      if(pFile == NULL)
      {
         LogError(ERR_FILE_IO, "Write(): Couldn't open OstModel.txt file");
         ExitProgram(1);
      }
      if(m_bMultiObjProblem == false)
      {
         fprintf(pFile,"Run   obj.function   ");
      }
      else
      {
         fprintf(pFile,"Run   ");
      }
      if(m_pObsGroup != NULL) m_pObsGroup->Write(pFile, WRITE_BNR, NULL);
      if(pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_BNR);
      m_pParamGroup->Write(pFile, WRITE_BNR);
      fprintf(pFile,"\n");
      fclose(pFile);
   }

   pFile = fopen(name, "a+");
	fprintf(pFile, "%-4d  ", m_Counter);

   if(m_bMultiObjProblem == false)
   {
      WritePreciseNumber(pFile, objFuncVal);
      fprintf(pFile, "  ");
   }

   if(m_pObsGroup != NULL) m_pObsGroup->Write(pFile, WRITE_SCI, m_CurMultiObjF);
   if(pRespVarGroup != NULL) pRespVarGroup->Write(pFile, WRITE_SCI);
   m_pParamGroup->Write(pFile, WRITE_SCI);
   fprintf(pFile, "\n");
   fclose(pFile);
} /* end Write() */

/******************************************************************************
WriteMetrics()
******************************************************************************/
void Model::WriteMetrics(FILE * pFile)
{  
   if(m_bUseSurrogates == true)
   {
      m_pDecision->WriteMetrics(pFile);
   }
   else
   {
      fprintf(pFile, "Total Evals             : %d\n", m_Counter);
      fprintf(pFile, "Telescoping Strategy    : ");
      switch(m_Telescope)
      {
         case(TSCOPE_PVEX): fprintf(pFile, "convex-power\n"); break;
         case(TSCOPE_CVEX): fprintf(pFile, "convex\n"); break;
         case(TSCOPE_LINR): fprintf(pFile, "linear\n"); break;
         case(TSCOPE_CAVE): fprintf(pFile, "concave\n"); break;
         case(TSCOPE_DCVE): fprintf(pFile, "delayed-concave\n"); break;
         default: fprintf(pFile, "none\n"); break;
      }
      if(m_bCaching == true)
         fprintf(pFile, "Cache Hits              : %d\n", m_NumCacheHits);
      if(m_pParameterCorrection != NULL)
         m_pParameterCorrection->WriteMetrics(pFile);
   }
} /* end WriteMetrics() */

/******************************************************************************
 * ExcludeConstantParameters()
 *
 * Checks that each parameter can take on multiple values. If not, a warning will
 * be reported and the parameter will be ignored by the optimizer..
 * ******************************************************************************/
void Model::ExcludeConstantParameters(void)
{
   if(m_pParamGroup == NULL) return;

   int i, j, nprm, id;
   double upr, lwr;
   char tmp1[DEF_STR_SZ];
   UnchangeableString prm_name;
   FILE * pLog;
   
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   if(id == 0)
   {
      pLog = fopen("OstConstants0.txt", "w");
   }
   else
   {
      pLog = NULL;
   }

   nprm = m_pParamGroup->GetNumParams();

   for(j = 0; j < nprm; j++)
   {
      prm_name = m_pParamGroup->GetParamPtr(j)->GetName();
      upr = m_pParamGroup->GetParamPtr(j)->GetUprBnd();
      lwr = m_pParamGroup->GetParamPtr(j)->GetLwrBnd();

      if(upr == lwr)
      {
         if(pLog != NULL)
         {
            fprintf(pLog, "Parameter %s will be treated as a constant\n", prm_name);
         }
         m_pParamGroup->ExcludeParam(prm_name);
         // ExcludeParam() modifies m_pParamGroup so restart the iterator
         j=0;
         nprm = m_pParamGroup->GetNumParams();
      }
   }
   if(pLog != NULL)
   {
      fclose(pLog);
   }
} /* end ExcludeConstantParameters() */

/******************************************************************************
CheckGlobalSensitivity()

Checks that each observation is sensitive to at least one parameter over the
range of possible parameter values. If an observation is not sensitive to any
parameters, a warning will be reported and the observation will be ignored in 
the calibration.

Also checks the sensitivity of each parameter. If a given parameter does not
affect the objective function over the entire parameter range, a warning will
be reported and the parameter will be ignored in the calibration.
******************************************************************************/
void Model::CheckGlobalSensitivity(void)
{
   if((m_pObsGroup == NULL) || (m_pParamGroup == NULL)) return;
   if(m_bCheckGlobalSens == false) return;

   int i, j, nobs, nprm; 
   double upr, lwr, Fupr, Flwr;
   double * pInit, * obs_sum, * prm_sum;
   double * ObsUpr, * ObsLwr;
   char tmp1[DEF_STR_SZ];
   UnchangeableString * obs_names, * prm_names;

   nobs = m_pObsGroup->GetNumObs();
   obs_sum = new double[nobs];
   ObsUpr = new double[nobs];
   ObsLwr = new double[nobs];
   obs_names = new UnchangeableString[nobs];

   nprm = m_pParamGroup->GetNumParams();
   pInit = new double[nprm];
   prm_sum = new double[nprm];
   prm_names = new UnchangeableString[nprm];

   //save intitial parameters
   m_pParamGroup->ReadParams(pInit);
   
   for(i = 0; i < nobs; i++)
   {
      obs_sum[i] = 0.00;
      obs_names[i] = m_pObsGroup->GetObsPtr(i)->GetName();
   }
   
   for(j = 0; j < nprm; j++)
   {
      prm_names[j] = m_pParamGroup->GetParamPtr(j)->GetName();
      upr = m_pParamGroup->GetParamPtr(j)->GetUprBnd();
      m_pParamGroup->GetParamPtr(j)->SetEstVal(upr);
      Fupr = Execute();

      //store computed observation values
      for(i = 0; i < nobs; i++)
      {
         ObsUpr[i] = m_pObsGroup->GetObsPtr(i)->GetComputedVal(true, true);
      }

      lwr = m_pParamGroup->GetParamPtr(j)->GetLwrBnd();
      m_pParamGroup->GetParamPtr(j)->SetEstVal(lwr);
      Flwr = Execute();

      //store computed observation values
      for(i = 0; i < nobs; i++)
      {
         ObsLwr[i] = m_pObsGroup->GetObsPtr(i)->GetComputedVal(true, true);
      }

      //store change in objective function
      prm_sum[j] = fabs(Fupr - Flwr);

      //accumulate changes in observation values
      for(i = 0; i < nobs; i++)
      {
         obs_sum[i] += fabs(ObsUpr[i] - ObsLwr[i]);
      }

      //restore parameter values
      m_pParamGroup->WriteParams(pInit);
   }/* end for() */

   //perform parameter sensitivity checks
   for(j = 0; j < nprm; j++)
   {
      if(prm_sum[j] <= NEARLY_ZERO)
      {
         sprintf(tmp1, "%s appears to be insensitive and has been set to a constant value", prm_names[j]);
         LogError(ERR_INS_PARM, tmp1);
         m_pParamGroup->ExcludeParam(prm_names[j]);
      }
   }

   //perform observation sensitivity checks: phase 1, log the warning messages
   for(i = 0; i < nobs; i++)
   {
      if(obs_sum[i] <= NEARLY_ZERO)
      {
         sprintf(tmp1, "%s appears to be insensitive and has been excluded from the calibration", obs_names[i]);
         LogError(ERR_INS_OBS, tmp1);
         m_pObsGroup->ExcludeObs(obs_names[i]);
      }
   }

   delete [] pInit;
   delete [] obs_sum;
   delete [] prm_sum;
   delete [] ObsUpr;
   delete [] ObsLwr;
   delete [] prm_names;
   delete [] obs_names;
} /* end CheckGlobalSensitivity() */

