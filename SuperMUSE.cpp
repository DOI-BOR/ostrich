/******************************************************************************
File     : SuperMUSE.cpp
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

The SuperMUSE class interfaces Ostrich with the EPA SuperMUSE cluster computing 
system; namely via the Java-based RepeatTasker program.  The two programs 
interface via a simple file-based communication strategy, where the 
presence/absence of certain files serve as a handshaking mechanism 
(see SuperMUSE.h for details).

Version History
07-13-07    lsm   added copyright information and initial comments.
******************************************************************************/
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "SuperMUSE.h"
#include "Model.h"
#include "ParameterGroup.h"

#include "SuperMuseUtility.h"
#include "Utility.h"
#include "Exception.h"

//definition of sleep() is platform-dependent
#ifdef WIN32
   #include <Windows.h> //needed for Sleep() and GetWindowsDirectory() commands
#else
   #include <unistd.h> //needed for sleep() command
#endif

/******************************************************************************
CTOR

Constructs a SuperMUSE class using the given input file.
******************************************************************************/
SuperMUSE::SuperMUSE(FILE * pFile, ModelABC * pModel)
{
   char * line;
   char tmp[DEF_STR_SZ];
   const char * pFileName = "OstIn.txt";
   m_pModel = pModel;

   //set reasonable defaults
   strcpy(m_Server, "0101Prog");
   strcpy(m_TaskFile, "SMuseTaskFile.txt");
   strcpy(m_TempFile, "SMuseTempFile.txt");
   strcpy(m_SuccessFile, "SMuseSuccessFile.txt");
   strcpy(m_ErrorFile, "SMuseErrorFile.txt");
   strcpy(m_ScriptFile, "SMuseScriptFile.txt");
   strcpy(m_ArgsFile, "SMuseArgumentsFile.txt");
   strcpy(m_ClientDir, "Simulations");
   strcpy(m_ServerDir, "FRAMESv2/Simulations");
   m_MaxJobTime     = 120; // 2 hours
   m_TaskID = 0;

   /* -----------------------------------------------
   Parse the SuperMUSE section of the input file.
   ----------------------------------------------- */
   //check for critical entries
   FindToken(pFile, "BeginSuperMUSE", pFileName);
   FindToken(pFile, "EndSuperMUSE", pFileName);
   rewind(pFile);

   FindToken(pFile, "BeginSuperMUSE", pFileName);
   line = GetNxtDataLine(pFile, pFileName);
      
   while(strncmp(line, "EndSuperMUSE", 12) != 0)
   {
      /*
      Technically, AllocatorServer is not correct...
      The server variable refers to the Ostrich Tasker
      which may or may not be on the same host as the 
      CPU Allocator.

      However, keeping AllocatorServer for backwards
      compatibility.
      */
      if(strncmp(line, "AllocatorServer", 15) == 0)
      {
         sscanf(line, "%s %s", tmp, m_Server);
      }
      else if(strncmp(line, "OstrichTaskerHostName", 21) == 0)
      {
         sscanf(line, "%s %s", tmp, m_Server);
      }
      else if(strncmp(line, "TaskFile", 8) == 0)
      {
         sscanf(line, "%s %s", tmp, m_TaskFile);
      }
      else if(strncmp(line, "TempFile", 8) == 0)
      {
         sscanf(line, "%s %s", tmp, m_TempFile);
      }
      else if(strncmp(line, "SuccessFile", 11) == 0)
      {
         sscanf(line, "%s %s", tmp, m_SuccessFile);
      }
      else if(strncmp(line, "ErrorFile", 9) == 0)
      {
         sscanf(line, "%s %s", tmp, m_ErrorFile);
      }
      else if(strncmp(line, "ScriptFile", 10) == 0)
      {
         sscanf(line, "%s %s", tmp, m_ScriptFile);
      }
      else if(strncmp(line, "ArgumentsFile", 13) == 0)
      {
         sscanf(line, "%s %s", tmp, m_ArgsFile);
      }
      else if(strncmp(line, "ClientDir", 9) == 0)
      {
         sscanf(line, "%s %s", tmp, m_ClientDir);
      }
      else if(strncmp(line, "ServerDir", 9) == 0)
      {
         sscanf(line, "%s %s", tmp, m_ServerDir);
      }
      else if(strncmp(line, "MaxJobTime", 10) == 0)
      {
         sscanf(line, "%s %d", tmp, &m_MaxJobTime);
      }
      else
      {
         sprintf(tmp, "SuperMUSE(): unknown token |%s|", line);
         LogError(ERR_FILE_IO, tmp);
      }

      line = GetNxtDataLine(pFile, pFileName);
   }/* end while() */   
                      
   IncCtorCount();
} /* end SuperMUSE::CTOR */

/******************************************************************************
LoadEnvVars()

Reads a list of environment variables from the specified file.
******************************************************************************/
void SuperMUSE::LoadEnvVars(char * pEnvVarFile)
{
   int j;
   EnvVarList * pCur, * pEnd;
   char var[1000], line[1000], * pTok;
   FILE * pIn;

   pCur = NULL;
   pEnd = NULL;
   pTok = NULL;
   pIn = fopen(pEnvVarFile, "r");
   if(pIn == NULL)
   {
	   LogError(ERR_FILE_IO, "Couldn't open SuperMUSE environment variable file");
	   ExitProgram(1);
   }
   while(!feof(pIn))
   {
      fgets(line, 1000, pIn);
      MyTrim(line);
      MyStrLwr(line);
      //two ways that a set command could be formatted
      if((strncmp(line, "set ", 4) == 0) || (strncmp(line, "call set ", 9) == 0))
      {
         pCur = new EnvVarList;
         pCur->pNxt = NULL;
         pTok = strstr(line, "set ") + 4;
         j = ExtractColString(pTok, var, '=');
         strcpy(pCur->pVar, "%");
         strcat(pCur->pVar, var);
         strcat(pCur->pVar, "%");
         pTok += j;
         strcpy(pCur->pVal, pTok);

         if(m_pEnvVars == NULL)
         {
            m_pEnvVars = pCur;
            pEnd = pCur;
         }
         else //tack on to end of list
         {
            pEnd->pNxt = pCur;
            pEnd = pCur;
         }
      }/*end if() */
   }/* end while() */
   fclose(pIn);
} /* end LoadEnvVars() */

/******************************************************************************
UnloadEnvVars()

Destroys a list of environment variables.
******************************************************************************/
void SuperMUSE::UnloadEnvVars(void)
{
   EnvVarList * pCur, * pNxt;
   //JB Fixed bad null assignment occurring at end of list due to final condition check
   //for(pCur = m_pEnvVars; pCur != NULL; pCur=pCur->pNxt)
   for(pCur = m_pEnvVars; pCur != NULL;)
   {
      pNxt = pCur->pNxt;
      delete pCur;
      pCur = pNxt;
   }
} /* end UnloadEnvVars() */

/******************************************************************************
ReplaceEnvVars()

Destroys a SuperMUSE class.
******************************************************************************/
void SuperMUSE::ReplaceEnvVars(char * pPathStr)
{
   EnvVarList * pCur;
   MyStrLwr(pPathStr);
   for(pCur = m_pEnvVars; pCur != NULL; pCur=pCur->pNxt)
      MyStrRep(pPathStr, pCur->pVar, pCur->pVal);
} /* end SuperMUSE::ReplaceEnvVars() */

/******************************************************************************
DTOR

Destroys a SuperMUSE class.
******************************************************************************/
void SuperMUSE::Destroy(void)
{   
   IncDtorCount();
} /* end SuperMUSE::DTOR */

/******************************************************************************
WriteTask()

Create an entry in the task file. This will list the script file name, plus
a set of command line arguments:
arg1 - server name
arg2 - task identifier
arg3 - file containing list of tasks
arg4 - client-side working directory
arg5 - server-side working directory
******************************************************************************/
void SuperMUSE::WriteTask(ParameterGroup * pGroup)
{
   FILE * pArgs;
   FILE * pTasks;

   pTasks = fopen(m_TempFile, "a");
   fprintf(pTasks, "%s %s %d %s %s %s \n", m_ScriptFile, m_Server, m_TaskID, 
                                          m_ArgsFile, m_ClientDir, m_ServerDir);
   fclose(pTasks);

   pArgs  = fopen(m_ArgsFile, "a");
   pGroup->WriteSuperMuseArgs(pArgs);
   fclose(pArgs);

   m_TaskID++;
}/* end WriteTask() */

/******************************************************************************
FinishTaskFile()

Finish up writing the temporary version of the task file, delete the error and
success files and rename the temporary task file to the actual file name 
required by the RepeatTasker program. This will launch the desired job on 
SuperMUSE in a safe manner.
******************************************************************************/
void SuperMUSE::FinishTaskFile(void)
{
   FILE * pFile;
   pFile = fopen(m_TempFile, "a");
   fprintf(pFile, "end\n");
   fclose(pFile);
   
   remove(m_SuccessFile);
   remove(m_ErrorFile);
   rename(m_TempFile, m_TaskFile);
}/* end FinishTaskFile() */

/******************************************************************************
WaitForTasker()

Wait for SuperMUSE to report back (via the success or error files).
******************************************************************************/
bool SuperMUSE::WaitForTasker(void)
{
   FILE * pSuccess;
   FILE * pError;
   int nticks;

   nticks = 0;
   while(1)
   {
#ifdef WIN32
      Sleep(1000); //1 second polling intervals
#else
      sleep(1); //1 second polling intervals
#endif

      nticks++;

      //check for successful completion
      pSuccess = fopen(m_SuccessFile, "r");
      if(pSuccess != NULL)
      {
         fclose(pSuccess);

         /* -------------------------
         Prepare for next go-around
         by deleting the args, task 
         and temp files and resetting
         the TaskID
         ------------------------- */
         remove(m_ArgsFile);
         remove(m_TempFile);
         remove(m_TaskFile);
         m_TaskID = 0;

         return true;
      }/* end if() */

      //check for failed completion
      pError = fopen(m_ErrorFile, "r");
      if(pError != NULL)
      {
         fclose(pError);
         LogError(ERR_SMUSE, "SuperMUSE Tasker failed to complete one or more tasks.");
         return false;            
      }/* end if() */

      if(nticks > (m_MaxJobTime*60))
      {
         LogError(ERR_SMUSE, "Timed out waiting for SuperMUSE Tasker.");
         return false;
      }/* end if() */

      //check for user abort
      if(IsQuit())
      {
         LogError(ERR_SMUSE, "User aborted SuperMUSE operation.");
         return false;
      }
   }/* end while() */         
 }/* end WaitForTasker() */

/******************************************************************************
GatherResult()

Collect the result (objective function value) of a task successfully 
completed on SuperMUSE. As each client finishes a task, it should create a 
directory named "TaskN" on the host compuer and directory under which Ostrich 
is running, where 'N' is the task number (defined as an input to the clients 
scriptfile). Each client is expected to copy model output files back to this
'TaskN' folder before signaling task completion. The GatherResult() routine
reads the model output files stored in the 'Nth' directory and computes and
returns an objective function for the associated model evaluation.
******************************************************************************/
double SuperMUSE::GatherResult(int taskid)
{
   double val;
   char taskDir[80];

   /* ----------------------------------
   Check to see if user has requested 
   program termination. If so, stuff
   return values with very large number
   and return.
   ---------------------------------- */
   if(IsQuit() == true)
   {
      return NEARLY_HUGE;
   }

   //gather result from the desired task   
   sprintf(taskDir, "Task%d", taskid);
   val = ((Model *)m_pModel)->GatherTask(taskDir);
   return val;
}/* end GatherResult() */

/******************************************************************************
WriteSetup()

Write SuperMUSE configuration to file.
******************************************************************************/
void SuperMUSE::WriteSetup(FILE * pFile)
{
   fprintf(pFile, "SuperMUSE Setup\n");
   fprintf(pFile, "SuperMUSE is ");
   if(IsSuperMUSE() == true) fprintf(pFile, "enabled\n");
   else                      fprintf(pFile, "disabled\n");
   fprintf(pFile, "Ostrich Tasker Host    : %s\n", m_Server);
   fprintf(pFile, "Task File              : %s\n", m_TaskFile);
   fprintf(pFile, "Temp File              : %s\n", m_TempFile);
   fprintf(pFile, "Success File           : %s\n", m_SuccessFile);
   fprintf(pFile, "Error File             : %s\n", m_ErrorFile);
   fprintf(pFile, "Script File            : %s\n", m_ScriptFile);
   fprintf(pFile, "Arguments File         : %s\n", m_ArgsFile);
   fprintf(pFile, "Client Working Folder  : %s\n", m_ClientDir);
   fprintf(pFile, "Server Working Folder  : %s\n", m_ServerDir);
   fprintf(pFile, "Max job Time (minutes) : %d\n\n", m_MaxJobTime);
}/* end WriteSetup() */

/******************************************************************************
EnvVarCleanup()

Cleanup environment variables.
******************************************************************************/
void SuperMUSE::EnvVarCleanup(void)
{
   //replace environment vars used in temp, task, success, error files
   m_pEnvVars = NULL;
   char iemVarsFile[1000];
   #ifdef WIN32
   GetWindowsDirectoryA(iemVarsFile, 1000);
   #else
   strcpy(iemVarsFile, "");
   #endif
   strcat(iemVarsFile, "\\iemSetCmdEnvironment.bat");
   LoadEnvVars(iemVarsFile); //create list of env vars
   ReplaceEnvVars(m_TempFile);
   ReplaceEnvVars(m_TaskFile);
   ReplaceEnvVars(m_SuccessFile);
   ReplaceEnvVars(m_ErrorFile);
   UnloadEnvVars(); //clean up list of env vars

   //clear out files, if they exist   
   remove(m_ArgsFile);
   remove(m_TempFile);
   remove(m_TaskFile);
} /* EnvVarCleanup() */
