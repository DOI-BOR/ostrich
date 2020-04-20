/******************************************************************************
File     : MemoryTracker.cpp
Author   : L. Shawn Matott
Copyright: 2013, L. Shawn Matott

Track memory usage of sub-routines.

Version History
01-29-13    lsm   added copyright information and initial comments.
******************************************************************************/
#include "MemoryTracker.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
void LogMemUsage(char * file, const char * tag)
{
	return;
}

double GetMemUsage(void)
{
	return -1.00;
}
#else

#ifdef WIN32
  #include "windows.h"
  #include "psapi.h"
#else
  #include "sys/types.h"
  #include "sys/sysinfo.h"
#endif

extern "C" {   
#ifndef WIN32
   double parseLine(char* line);
   double getVmSize(void);
   double getVmRSS(void);
#endif
}

/******************************************************************************
GetMemUsage()

Retrieve phyiscal memory usage.

Code for Windows and Linux was taken from the following site: 
http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
******************************************************************************/
double GetMemUsage(void)
{  
#ifdef WIN32
   MEMORYSTATUSEX memInfo;
   memInfo.dwLength = sizeof(MEMORYSTATUSEX);
   GlobalMemoryStatusEx(&memInfo);
   double totalVirtualMem = (double)(memInfo.ullTotalPageFile);
   double virtualMemUsed = (double)(memInfo.ullTotalPageFile - memInfo.ullAvailPageFile);

   PROCESS_MEMORY_COUNTERS pmc;
   GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
   double virtualMemUsedByMe = (double)(pmc.PagefileUsage);

   double totalPhysMem = (double)(memInfo.ullTotalPhys);

   double physMemUsed = (double)(memInfo.ullTotalPhys - memInfo.ullAvailPhys);

   double physMemUsedByMe = (double)(pmc.WorkingSetSize);
#elif __APPLE__
	double physMemUsedByMe = -1.00;
#else
   struct sysinfo memInfo;
    
   sysinfo (&memInfo);
   double totalVirtualMem = memInfo.totalram;
   //Add other values in next statement to avoid int overflow on right hand side...
   totalVirtualMem += memInfo.totalswap;
   totalVirtualMem *= memInfo.mem_unit;

   double virtualMemUsed = memInfo.totalram - memInfo.freeram;
   //Add other values in next statement to avoid int overflow on right hand side...
   virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
   virtualMemUsed *= memInfo.mem_unit;

   double virtualMemUsedByMe = (double)(getVmSize()*1000);

   double totalPhysMem = memInfo.totalram;
   //Multiply in next statement to avoid int overflow on right hand side...
   totalPhysMem *= memInfo.mem_unit;

   double physMemUsed = memInfo.totalram - memInfo.freeram;
   //Multiply in next statement to avoid int overflow on right hand side...
   physMemUsed *= memInfo.mem_unit;

   double physMemUsedByMe = (double)(getVmRSS()*1000);
#endif
   return physMemUsedByMe;
}/* end GetMemUsage() */

/******************************************************************************
LogMemUsage()

Log the current memory usage to a file. If 'file' argument is NULL print to 
stdout instead.

Code for Windows and Linux was taken from the following site: 
http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
******************************************************************************/
void LogMemUsage(char * file, const char * tag)
{  
   FILE * pFile;   
   if(file == NULL) pFile = stdout;
   else pFile = fopen(file, "a");

#ifdef WIN32
   MEMORYSTATUSEX memInfo;
   memInfo.dwLength = sizeof(MEMORYSTATUSEX);
   GlobalMemoryStatusEx(&memInfo);
   double totalVirtualMem = (double)(memInfo.ullTotalPageFile);
   fprintf(pFile, "%s,TotalVirtualMemory,%.0lf\n", tag, totalVirtualMem);

   double virtualMemUsed = (double)(memInfo.ullTotalPageFile - memInfo.ullAvailPageFile);
   fprintf(pFile, "%s,VirtualMemoryInUse,%.0lf\n", tag, virtualMemUsed);

   PROCESS_MEMORY_COUNTERS pmc;
   GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
   double virtualMemUsedByMe = (double)(pmc.PagefileUsage);
   fprintf(pFile, "%s,VirtualMemoryInUseByMe,%.0lf\n", tag, virtualMemUsedByMe);

   double totalPhysMem = (double)(memInfo.ullTotalPhys);
   fprintf(pFile, "%s,TotalPhysicalMemory,%.0lf\n", tag, totalPhysMem);

   double physMemUsed = (double)(memInfo.ullTotalPhys - memInfo.ullAvailPhys);
   fprintf(pFile, "%s,PhysicalMemoryInUse,%.0lf\n", tag, physMemUsed);

   double physMemUsedByMe = (double)(pmc.WorkingSetSize);
   fprintf(pFile, "%s,PhysicalMemoryInUseByMe,%.0lf\n", tag, physMemUsedByMe);
#else
   struct sysinfo memInfo;
    
   sysinfo (&memInfo);
   double totalVirtualMem = memInfo.totalram;
   //Add other values in next statement to avoid int overflow on right hand side...
   totalVirtualMem += memInfo.totalswap;
   totalVirtualMem *= memInfo.mem_unit;
   fprintf(pFile, "%s,TotalVirtualMemory,%.0lf\n", tag, totalVirtualMem);

   double virtualMemUsed = memInfo.totalram - memInfo.freeram;
   //Add other values in next statement to avoid int overflow on right hand side...
   virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
   virtualMemUsed *= memInfo.mem_unit;
   fprintf(pFile, "%s,VirtualMemoryInUse,%.0lf\n", tag, virtualMemUsed);

   double virtualMemUsedByMe = (double)getVmSize()*1000;
   fprintf(pFile, "%s,VirtualMemoryInUseByMe,%.0lf\n", tag, virtualMemUsedByMe);

   double totalPhysMem = memInfo.totalram;
   //Multiply in next statement to avoid int overflow on right hand side...
   totalPhysMem *= memInfo.mem_unit;
   fprintf(pFile, "%s,TotalPhysicalMemory,%.0lf\n", tag, totalPhysMem);

   double physMemUsed = memInfo.totalram - memInfo.freeram;
   //Multiply in next statement to avoid int overflow on right hand side...
   physMemUsed *= memInfo.mem_unit;
   fprintf(pFile, "%s,PhysicalMemoryInUse,%.0lf\n", tag, physMemUsed);

   double physMemUsedByMe = (double)getVmRSS()*1000;
   fprintf(pFile, "%s,PhysicalMemoryInUseByMe,%.0lf\n", tag, physMemUsedByMe);
#endif
   if(file != NULL) fclose(pFile);
}/* end ShowPhysMemUsage() */


#ifndef WIN32
   double parseLine(char* line)
   {
      int i = strlen(line);
      while (*line < '0' || *line > '9') line++;
      line[i-3] = '\0';
      double f = atof(line);
      return f;
   }

   //Note: this value is in KB!
   double getVmSize(void)
   { 
      FILE* file = fopen("/proc/self/status", "r");
      int result = -1;
      char line[128];
    
      while (fgets(line, 128, file) != NULL)
      {
          if (strncmp(line, "VmSize:", 7) == 0)
          {
              result = parseLine(line);
              break;
          }
      }
      fclose(file);
      return result;
    }

    //Note: this value is in KB!
    double getVmRSS(void)
    { 
       FILE* file = fopen("/proc/self/status", "r");
       int result = -1;
       char line[128];
    
       while (fgets(line, 128, file) != NULL)
       {
           if (strncmp(line, "VmRSS:", 6) == 0)
           {
               result = parseLine(line);
               break;
           }
       }
       fclose(file);
       return result;
    }
#endif
#endif
