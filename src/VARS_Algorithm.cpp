/******************************************************************************
File       : VARS_Algorithm.cpp
Author     : L. Shawn Matott
Copyrights : 2018, L. Shawn Matott 

VARS - Variogram Analysis of Response Surfaces 

Interface class for the VARS plugin, which must be licensed separately.
 
Version History
06-18-18    lsm   created file
******************************************************************************/
#include <mpi.h>
#include <string.h>
#include "VARS_Algorithm.h"
#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"
#include "Model.h"
#include "ParameterGroup.h"
#include "ObjectiveFunction.h"

#ifdef WIN32
  #include <Windows.h>
#else
  #include <dlfcn.h>
#endif

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
VARS::VARS(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_PluginFile = NULL;
   m_Subroutine = NULL;
   
   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by PAES and it's member variables.
******************************************************************************/
void VARS::Destroy(void)
{
   delete [] m_PluginFile;
   delete [] m_Subroutine;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

In the VARS contextm Calibrate() is synonymous with Optimize()
******************************************************************************/
void VARS::Calibrate(void)
{
   Optimize();
} /* end Calibrate() */

/******************************************************************************
VARS::InitFromFile()

Reads configurations parameters from a file
******************************************************************************/
void VARS::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];

   //read in VARS configuration
   pFile = fopen(pFileName, "r");
   if (pFile == NULL)
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open VARS config. file.");
      ExitProgram(-1);
      return;
   }/* end if() */ 

   //make sure correct tokens are present
   if (CheckToken(pFile, "Begin_VARS_Plugin", pFileName) == true)
   {
      FindToken(pFile, "End_VARS_Plugin", pFileName);
      rewind(pFile);

      FindToken(pFile, "Begin_VARS_Plugin", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while (strstr(line, "End_VARS_Plugin") == NULL)
      {
         if (strstr(line, "LibraryFile") != NULL)
         {
            m_PluginFile = new char[strlen(line)];
            strcpy(m_PluginFile, &line[11]);
            MyTrim(m_PluginFile);
         }/*end else if() */
         else if (strstr(line, "Subroutine") != NULL)
         {
            m_Subroutine = new char[strlen(line)];
            strcpy(m_Subroutine, &line[10]);
            MyTrim(m_Subroutine);
         }/*end else if() */
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void VARS::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Variogram Analysis of Response Surfaces\n");
   fprintf(pFile, "VARS Library File       : %s\n", m_PluginFile);
   fprintf(pFile, "VARS Subroutine         : %s\n", m_Subroutine);
  
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
VARS::Optimize()

Run the selected plugin routine.
******************************************************************************/
void VARS::Optimize(void)
{
   #ifdef WIN32
      HINSTANCE lib_handle;
   #else  
      void * lib_handle = NULL;
      char * error;
   #endif

   /* plugin function pointers */
   int (*BoxCoxDevTest)(int, char **);

   char msg[DEF_STR_SZ];
   StatusStruct pStatus;

   InitFromFile(GetInFileName());

   WriteSetup(m_pModel, "VARS - Variogram Analysis of Response Surfaces");
   //write banner
   WriteBanner(m_pModel, "samp  ", "Precent Complete");

   /* attempt to import shared library or DLL */
   #ifdef WIN32
      lib_handle = LoadLibraryA(m_PluginFile);

      if(!lib_handle)
      {
         sprintf(msg, "VARS_Algorithm :: LoadLibrary(%s) failed", m_PluginFile);
         LogError(ERR_BAD_ARGS, msg);
         ExitProgram(-1);
      }
   #else
      lib_handle = dlopen(m_PluginFile, RTLD_LAZY);
      if(!lib_handle)
      {
         sprintf(msg, "VARS_Algorithm :: dlopen(%s) failed with error (%s)", m_PluginFile, dlerror());
         LogError(ERR_BAD_ARGS, msg);
         ExitProgram(-1);
      }
   #endif

   /* load the symbols */
   #ifdef WIN32
      BoxCoxDevTest = (int(*)(int, char**))GetProcAddress(lib_handle, "VARS_main");
      if(BoxCoxDevTest == NULL)
      {
         sprintf(msg,"VARS_Algorithm :: failed to load DLL function (VARS_main)");
         LogError(ERR_BAD_ARGS, msg);
         ExitProgram(-1);
      }        
   #else
      dlerror();
      BoxCoxDevTest = (int(*)(int, char**))dlsym(lib_handle, "VARS_main");
      error = dlerror();
      if(error != NULL)
      {
         sprintf(msg,"VARS_Algorithm :: failed to load DLL function (VARS_main) - error code: %s", error);
         LogError(ERR_BAD_ARGS, msg);
         ExitProgram(-1);
      }
   #endif

   /* which subroutine? */
   if(m_Subroutine == NULL)
   {
      sprintf(msg, "No VARS subroutine was specified");
      LogError(ERR_BAD_ARGS, msg);
      ExitProgram(-1);
   }
   else if(strcmp(m_Subroutine, "BoxCoxDevTest") == 0)
   {
      int x = (*BoxCoxDevTest)(0, NULL);
   }
   else if(strcmp(m_Subroutine, "main_VARS") == 0)
   {
   }
   else
   {
      sprintf(msg, "Unkown VARS subroutine (%s)", m_Subroutine);
      LogError(ERR_BAD_ARGS, msg);
      ExitProgram(-1);
   }

   /* final status */
   pStatus.curIter = m_CurIter;
   pStatus.maxIter = m_CurIter;
   pStatus.pct = (float)100.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   //write algorithm metrics
   WriteAlgMetrics(this);

   /* unload the library */
   #ifdef WIN32
      FreeLibrary(lib_handle);
   #else
      dlclose(lib_handle);
   #endif

} /* end Optimize() */

/******************************************************************************
VARS_Program()

Run the desired VARS subroutine.
******************************************************************************/
void VARS_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("VARS", 1);
   VARS * TheAlg = new VARS(model);
   MEM_CHECK(TheAlg);

   if (model->GetObjFuncId() == OBJ_FUNC_WSSE) { TheAlg->Calibrate(); }
   else { TheAlg->Optimize(); }

   delete TheAlg;
   delete model;
} /* end VARS_Program() */

