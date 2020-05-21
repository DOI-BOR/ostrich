/******************************************************************************
File      : file_mpi.c
Author    : L. Shawn Matott
Copyright : 2013, L. Shawn Matott

Implements mpi functions using file-based approach so that one can compile 
MPI code in environments that don't have MPI libraries.

Version History
10-24-13    lsm   added copyright information and initial comments.
******************************************************************************/
#include "file_mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#ifdef WIN32
  #include <windows.h>
  #include <direct.h>
#endif

#ifdef USE_FILE_MPI /* change preprocessor def if no file-based mpi desired */

/* global variables */
int FMPI_gMpiSize = 1;
int FMPI_gMpiRank = 0;
int FMPI_gMyPid = -1;
int FMPI_gMpiIsInitialized = 0;
unsigned int FMPI_gTimeout = 60000;
unsigned int ** FMPI_gMsgIds = NULL;

/* max string buffers */
#define FMPI_MAX_LINE_SIZE (1024)
#define FMPI_MAX_FNAME_SIZE (1024)

/* timeouts for initial processor synchronization */
#define FMPI_SYNC_TIMEOUT_MS (FMPI_gTimeout)
#define FMPI_SYNC_POLL_INTERVAL_MS (1000)

/* timeouts for MPI_Barrier() */
#define FMPI_BARRIER_TIMEOUT_MS (FMPI_gTimeout)
#define FMPI_BARRIER_POLL_INTERVAL_MS (10)

/* timeouts for MPI_Finalize() */
#define FMPI_FINALIZE_TIMEOUT_MS (FMPI_gTimeout)
#define FMPI_FINALIZE_POLL_INTERVAL_MS (10)

/* timeouts for opening files for reading */
#define FMPI_RECV_TIMEOUT_MS (FMPI_gTimeout)
#define FMPI_RECV_POLL_INTERVAL_MS (10)

/* timeouts for deleting files */
#define FMPI_DEL_TIMEOUT_MS (FMPI_gTimeout)
#define FMPI_DEL_POLL_INTERVAL_MS (10)

/* error code for aborting in response to aborts by other ranks */
#define MPI_ERR_ABORT_DETECTED (-999)

/* private helper functions */
int FMPI_CountFiles(char * prefix);
int FMPI_CountSendTags(void);
int FMPI_GetSendTagList(int * pList, int nTags);
int FMPI_GetProcessList(unsigned int * pList, int n, char * prefix);
int FMPI_GetPid(void);
void FMPI_Sleep(int millisecs);
void FMPI_RemoveFile(char * fname);
int FMPI_GetTypeSize(MPI_Datatype type);
int * FMPI_SortSources(unsigned int ** msgIds, int n, int dst);

/* private helper functions to synch send/recv operations */
FILE * FMPI_OpenFileForSend(char * fname);
void FMPI_CloseFileForSend(FILE * pFile, char * fname);
FILE * FMPI_OpenFileForRecv(char * fname);
int FMPI_FileExists(char * fname);

/* private helper functions for reductions */
int FMPI_IntegerReduction(int cur, int val, MPI_Op op);
double FMPI_DoubleReduction(double cur, double val, MPI_Op op);
char FMPI_CharaceterReduction(char cur, char val, MPI_Op op);

/* private helper that checks for Abort */
void FMPI_CheckForAbort(void);

/********************************************************************
FMPI_CheckForAbort()

Check for presense of MPI_Abort.txt and abort if needed.
********************************************************************/
void FMPI_CheckForAbort(void)
{
   FILE * pAbort = fopen("MPI_Abort.txt", "r");
   if(pAbort != NULL)
   {
      fclose(pAbort);
      MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ABORT_DETECTED);
   }
}/* end FMPI_CheckForAbort() */

/********************************************************************
FMPI_GetTypeSize()

Get number of bytes in datatype.
********************************************************************/
int FMPI_GetTypeSize(MPI_Datatype type)
{
   if(type == MPI_CHAR) return sizeof(char);
   if(type == MPI_INTEGER) return sizeof(int);
   if(type == MPI_DOUBLE) return sizeof(double);
   return 0;
}/* end FMPI_GetTypeSize() */

/********************************************************************
FMPI_RemoveFile()

Attempt to remove/delete a file.
********************************************************************/
void FMPI_RemoveFile(char * fname)
{
   FILE * pFile;
   unsigned int howLong = 0;
   while(howLong < FMPI_DEL_TIMEOUT_MS)
   {
      FMPI_CheckForAbort();

      remove(fname);

      pFile = fopen(fname, "r");
      if(pFile == NULL)
      {
         //printf("rank # %d deleted a file (%s)!\n", FMPI_gMpiRank, fname);
         return;
      }
      fclose(pFile);

      FMPI_Sleep(FMPI_DEL_POLL_INTERVAL_MS);
      howLong += FMPI_DEL_POLL_INTERVAL_MS;
   }/* end while(waiting at barrier) */

   printf("Error - rank # %d timed out trying to delete a file (%s)!\n", FMPI_gMpiRank, fname);
   MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
}/* end FMPI_RemoveFile() */

/********************************************************************
FMPI_IntegerReduction()

Perform a reduction operation on an integer.
********************************************************************/
int FMPI_IntegerReduction(int cur, int val, MPI_Op op)
{
   switch(op)
   {
      case(MPI_SUM) :
      {
         return(cur+val);
      }
      case(MPI_MIN) :
      {
         if(cur < val) return cur;
         return val;
      }
      case(MPI_MAX) :
      {
         if(cur > val) return cur;
         return val;
      }
      default :
      {
         return cur;
      }
   }/* end switch() */
   return cur;
}/* end FMPI_IntegerReduction() */

/********************************************************************
FMPI_DoubleReduction()

Perform a reduction operation on a double.
********************************************************************/
double FMPI_DoubleReduction(double cur, double val, MPI_Op op)
{
   switch(op)
   {
      case(MPI_SUM) :
      {
         return(cur+val);
      }
      case(MPI_MIN) :
      {
         if(cur < val) return cur;
         return val;
      }
      case(MPI_MAX) :
      {
         if(cur > val) return cur;
         return val;
      }
      default :
      {
         return cur;
      }
   }/* end switch() */
   return cur;
}/* end FMPI_DoubleReduction() */

/********************************************************************
FMPI_CharacterReduction()

Perform a reduction operation on a character.
********************************************************************/
char FMPI_CharacterReduction(char cur, char val, MPI_Op op)
{
   switch(op)
   {
      case(MPI_SUM) :
      {
         return(cur+val);
      }
      case(MPI_MIN) :
      {
         if(cur < val) return cur;
         return val;
      }
      case(MPI_MAX) :
      {
         if(cur > val) return cur;
         return val;
      }
      default :
      {
         return cur;
      }
   }/* end switch() */
   return cur;
}/* end FMPI_CharacterReduction() */

/********************************************************************
FMPI_Sleep()

Wait a bit.
********************************************************************/
void FMPI_Sleep(int millisecs)
{
   #ifdef WIN32
      Sleep(millisecs);
   #else
      sleep(millisecs);
   #endif
}/* end FMPI_Sleep() */

/********************************************************************
FMPI_GetPid()

Retrieve a unique process/thread id.
********************************************************************/
int FMPI_GetPid(void)
{
   //already assigned?
   if(FMPI_gMyPid != -1)
      return FMPI_gMyPid;

   #ifdef WIN32
      FMPI_gMyPid = GetCurrentThreadId();
   #else
     FMPI_gMyPid = getpid();
   #endif

   return FMPI_gMyPid;
}/* end FMPI_GetPid() */

/********************************************************************
FMPI_SortSources()

Sort processors in ascending order based on message id. Sort is 
performed using insertion sort. For details see: 
   http://en.wikipedia.org/wiki/Insertion_sort
********************************************************************/
int * FMPI_SortSources(unsigned int ** msgIds, int n, int dst)
{
   int i,s,e,j,t;
   int * out = (int *)(malloc(sizeof(int)*n));
   int * val = (int *)(malloc(sizeof(int)*n));

   //initialize
   for(i = 0; i < n; i++)
   {
      out[i] = i;
      val[i] = msgIds[i][dst];
   }

   //insertion sort
   for(s = 1; s < n; s++)
   {
      for(e = 0; e < s; e++)
      {
         t = val[s];
         if(t < val[e])
         {
            for(j = s; j > e; j--)
            {
               val[j] = val[j-1];
               out[j] = out[j-1];
            }
            val[e] = t;
            out[e] = s;
         }
      }
   }
   free(val);
   return out;
}/* end FMPI_SortSources() */

/********************************************************************
MPI_Wtime()

Retrieve elapsed time.
********************************************************************/
double MPI_Wtime(void)
{
   return((double)time(NULL));
}/* end MPI_Wtime() */

/********************************************************************
MPI_Get_processor_name()

Use hostname command to obtain processor name.
********************************************************************/
int MPI_Get_processor_name(char *name, int *resultlen)
{
   int pid = FMPI_GetPid();
   int i,j;
   char tmpFileName[FMPI_MAX_FNAME_SIZE];
   char cmdLine[FMPI_MAX_LINE_SIZE];
   char tmpStr[MPI_MAX_PROCESSOR_NAME];
   FILE * pOut;

   sprintf(tmpFileName, "runFileMPI.hostname.%d", pid);
   sprintf(cmdLine, "hostname > %s", tmpFileName);
   name[0] = (char)NULL;
   system(cmdLine);
   pOut = fopen(tmpFileName, "r");
   if(pOut == NULL) return MPI_ERROR;
   fgets(tmpStr, MPI_MAX_PROCESSOR_NAME, pOut);
   fclose(pOut);
   FMPI_RemoveFile(tmpFileName);

   j=0;
   for(i = 0; i < (int)strlen(tmpStr); i++)
   {
      if((tmpStr[i] != ' ')  && 
         (tmpStr[i] != '\r') && 
         (tmpStr[i] != '\n')) 
      {
         name[j++] = tmpStr[i];
      }
   }
   name[j] = (char)NULL;
   *resultlen = strlen(name);
   return MPI_SUCCESS;
}/* end MPI_Get_processor_name() */

/********************************************************************
MPI_Init()

Initialize rank and size by processing config file.
********************************************************************/
int MPI_Init(int * argc, char *** argv)
{
   int i, j, num;
   unsigned int * pProcIds;
   char syncFileName[FMPI_MAX_FNAME_SIZE];
   FILE * pSync;
   unsigned int mypid;
   unsigned int howLong = 0;
   int numSync = 0;
   char line[FMPI_MAX_LINE_SIZE];
   FILE * pCfg;

   if(FMPI_gMpiIsInitialized == 1)
   {
      printf("Error --- MPI already initialized!\n");
      MPI_Barrier(MPI_COMM_WORLD); //synchronize processors
      return MPI_ERROR;
   }

   pCfg = fopen("FileMpiIn.txt", "r");
   if(pCfg == NULL)
   {
      printf("Error --- could not open FileMPI configuration file (FileMpiIn.txt)\n");
      MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
      return MPI_ERROR;
   }

   //set confugration defaults
   FMPI_gTimeout = 30;
   FMPI_gMpiSize = 1;

   //parse configuration file
   while(!feof(pCfg))
   {
      fgets(line, FMPI_MAX_LINE_SIZE, pCfg);
      if(strncmp(line, "Timeout", strlen("Timeout")) == 0)
      {
         FMPI_gTimeout = atoi(&(line[strlen("Timeout")]));
      }
      else if(strncmp(line, "NumCores", strlen("NumCores")) == 0)
      {
         FMPI_gMpiSize = atoi(&(line[strlen("NumCores")]));
      }
   }/* end while() */
   fclose(pCfg);

   if(FMPI_gMpiSize < 1)
   {
      printf("Error --- invalid NumCores entry in FileMPI configuration file\n");
      printf("      --- (%s)\n", line);
      MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
      return MPI_ERROR;
   }/* end if() */

   /* ---------------------------------------------------------------------------
   Synchronize the various cores using files.
     First create a file containing the process id (runFileMPI.Sync.[pid])

     Then wait for all other processors to write their version of the file.

     Then read in the list of process id and assign ranks based on sorted 
     pids.

     Throw an error if not all processors create a sync file within some 
     timeout period.

     Delete sync files when finished.
   --------------------------------------------------------------------------- */
   mypid = FMPI_GetPid();
   sprintf(syncFileName, "runFileMPI.Sync.%d", mypid);
   pSync = fopen(syncFileName, "w");
   fprintf(pSync, "%d", mypid);
   fclose(pSync);

   while((numSync = FMPI_CountFiles("runFileMPI.Sync.")) < FMPI_gMpiSize)
   {
      FMPI_CheckForAbort();
      FMPI_Sleep(FMPI_SYNC_POLL_INTERVAL_MS);
      howLong += FMPI_SYNC_POLL_INTERVAL_MS;
      if(howLong >= FMPI_SYNC_TIMEOUT_MS)
      {
         printf("Error - timed out waiting for processors to sync up!\n");
         MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
         break;
      }/* end if() */
   }/* end while() */

   if(numSync < FMPI_gMpiSize)
   {
      printf("Only %d out of %d are communicating.\n", numSync, FMPI_gMpiSize);
      FMPI_RemoveFile(syncFileName);
      MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
      return MPI_ERROR;
   }/* end if (one or more processors not reporting) */
        
   /* assign MPI rank based on rank of process id */
   pProcIds = (unsigned int *)(malloc(FMPI_gMpiSize*sizeof(unsigned int)));
   num = FMPI_GetProcessList(pProcIds, FMPI_gMpiSize, "runFileMPI.Sync.");
   /* error check the actual list size */
   if(num != FMPI_gMpiSize)
   {
      printf("Error - size mismatch in GetProcessList() (%d vs. %d)\n", num, FMPI_gMpiSize);
      free(pProcIds);
      FMPI_RemoveFile(syncFileName);
      MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
      return MPI_ERROR;
   }/* end if() */
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      if(pProcIds[i] < mypid) FMPI_gMpiRank++;
   }/* end for() */
   //printf("Process ID = %d, MPI rank = %d\n", mypid, gMpiRank);
   free(pProcIds);

   MPI_Barrier(MPI_COMM_WORLD);
   FMPI_RemoveFile(syncFileName);
   FMPI_gMpiIsInitialized = 1;

   //cleanup configuration and launcher files
   remove("FileMpiRunWin32.bat");
   remove("FileMpiIn.txt");

   //create array of message ids
   FMPI_gMsgIds = (unsigned int **)(malloc(sizeof(unsigned int *)*FMPI_gMpiSize));
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      FMPI_gMsgIds[i] = (unsigned int *)(malloc(sizeof(unsigned int)*FMPI_gMpiSize));
   }
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      for(j = 0; j < FMPI_gMpiSize; j++)
      {
         FMPI_gMsgIds[i][j] = 0;
      }    
   }

   return MPI_SUCCESS;
}/* end MPI_Init() */

/********************************************************************
MPI_Abort()

Quit program.
********************************************************************/
int MPI_Abort(MPI_Comm comm, int errorcode)
{
   int i;
   FILE * pAbort;

   printf("Rank # %d has entered MPI_Abort() with an error code of %d\n", FMPI_gMpiRank, errorcode);
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      free(FMPI_gMsgIds[i]);
   }
   free(FMPI_gMsgIds);

   //create message file to let other procesors know something went wrong
   pAbort = fopen("MPI_Abort.txt", "a");
   fprintf(pAbort, "Rank # %d has aborted with an error code of %d!\n", FMPI_gMpiRank, errorcode);
   fclose(pAbort);

   exit(-1);
   return 0;
}/* end MPI_Abort() */

/********************************************************************
MPI_Comm_size()

Get number of processors.
********************************************************************/
int MPI_Comm_size(MPI_Comm comm, int * size)
{
   *size = FMPI_gMpiSize;
   return MPI_SUCCESS;
}/* end MPI_Comm_size() */

/********************************************************************
MPI_Comm_rank()

Get rank.
********************************************************************/
int MPI_Comm_rank(MPI_Comm comm, int * rank)
{
   *rank = FMPI_gMpiRank;
   return MPI_SUCCESS;
}/* end MPI_Comm_rank() */

/********************************************************************
MPI_Allgather()

"All gather" operation.
********************************************************************/
int MPI_Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               MPI_Comm comm)
{
   int i;
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, 
                 recvtype, i, comm);
   }/* end for() */
   return MPI_SUCCESS;
}/* end MPI_Allgather() */

/********************************************************************
MPI_Allgatherv()

Variable length and displacement "all gather" operation.
********************************************************************/
int MPI_Allgatherv (void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcounts, int *displs, 
                     MPI_Datatype recvtype, MPI_Comm comm )
{
   int i;
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, 
                  displs, recvtype, i, comm);
   }/* end for() */
   return MPI_SUCCESS;
}/* end MPI_Allgatherv() */

/********************************************************************
MPI_Gatherv()

Variable length and displacement gather operation.
********************************************************************/
int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                void *recvbuf, int *recvcnts, int *displs, 
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
   MPI_Status status;
   void * dst;
   int i;
   int offset;
   int typesize;
   int nbytes;

   typesize = FMPI_GetTypeSize(sendtype);
   
   if(FMPI_gMpiRank == root)
   {
      //copy root data from sendbuf to recvbuf
      offset = displs[root]*typesize;
      nbytes = sendcnt*typesize;
      dst = &(((char *)recvbuf)[offset]);
      memcpy(dst, sendbuf, nbytes);

      //execute sequence of receives from other ranks
      for(i = 0; i < FMPI_gMpiSize; i++)
      {
         if(i  != root)
         {
            offset = displs[i]*typesize;
            dst = &(((char *)recvbuf)[offset]);
            MPI_Recv(dst, recvcnts[i], recvtype, i, MPI_ANY_TAG, comm, &status);
         }/* end if() */
      }/* end for() */
   }/* end if() */
   else // send data to root
   {
      MPI_Send(sendbuf, sendcnt, sendtype, root, MPI_ANY_TAG, comm);
   }/* end if() */
   return MPI_SUCCESS;
}/* end MPI_Gatherv() */

/********************************************************************
MPI_Gather()

Gather operation.
********************************************************************/
int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm)
{
   MPI_Status status;
   void * dst;
   int i;
   int offset;
   int typesize;
   int nbytes;

   typesize = FMPI_GetTypeSize(sendtype);
   
   if(FMPI_gMpiRank == root)
   {
      //copy root data from sendbuf to recvbuf
      offset = root*sendcnt*typesize;
      nbytes = sendcnt*typesize;
      dst = &(((char *)recvbuf)[offset]);
      memcpy(dst, sendbuf, nbytes);

      //execute sequence of receives from other ranks
      for(i = 0; i < FMPI_gMpiSize; i++)
      {
         if(i  != root)
         {
            offset = i*sendcnt*typesize;
            dst = &(((char *)recvbuf)[offset]);
            MPI_Recv(dst, recvcnt, recvtype, i, MPI_ANY_TAG, comm, &status);
            //printf("Rank %d : Got gather message from %d\n", FMPI_gMpiRank, status.MPI_SOURCE);
         }/* end if() */
      }/* end for() */
   }/* end if() */
   else //send data to root
   {
      MPI_Send(sendbuf, sendcnt, sendtype, root, MPI_ANY_TAG, comm);
   }/* end if() */

   return MPI_SUCCESS;
}/* end MPI_Gather() */

/********************************************************************
MPI_Scatter()

Scatter operation.
********************************************************************/
int MPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
                int root, MPI_Comm comm)
{
   MPI_Status status;
   void * src;
   int i;
   int offset;
   int typesize;
   int nbytes;

   typesize = FMPI_GetTypeSize(sendtype);
   
   if(FMPI_gMpiRank == root)
   {
      //copy root data from sendbuf to recvbuf
      offset = root*recvcnt*typesize;
      nbytes = recvcnt*typesize;
      src = &(((char *)sendbuf)[offset]);
      memcpy(recvbuf, src, nbytes);

      //execute sequence of sends to other ranks
      for(i = 0; i < FMPI_gMpiSize; i++)
      {
         if(i  != root)
         {
            offset = i*recvcnt*typesize;
            src = &(((char *)sendbuf)[offset]);
            MPI_Send(src, sendcnt, sendtype, i, MPI_ANY_TAG, comm);
            //printf("Rank %d : Sent scatter message to %d\n", FMPI_gMpiRank, i);
         }/* end if() */
      }/* end for() */
   }/* end if() */
   else //receive data from root
   {
      MPI_Recv(recvbuf, recvcnt, recvtype, root, MPI_ANY_TAG, comm, &status);
   }/* end if() */

   return MPI_SUCCESS;
}/* end MPI_Scatter() */

/********************************************************************
MPI_Scatterv()

Variable length and displacement scatter operation.
********************************************************************/
int MPI_Scatterv(void *sendbuf, int * sendcnts, int *displs, 
                MPI_Datatype sendtype, void *recvbuf, int recvcnt, 
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
   MPI_Status status;
   void * src;
   int i;
   int offset;
   int typesize;
   int nbytes;

   typesize = FMPI_GetTypeSize(recvtype);
   
   if(FMPI_gMpiRank == root)
   {
      //copy root data from sendbuf to recvbuf
      offset = displs[root]*typesize;
      nbytes = recvcnt*typesize;
      src = &(((char *)sendbuf)[offset]);
      memcpy(recvbuf, src, nbytes);

      //execute sequence of sends to other ranks
      for(i = 0; i < FMPI_gMpiSize; i++)
      {
         if(i  != root)
         {
            offset = displs[i]*typesize;
            src = &(((char *)sendbuf)[offset]);
            MPI_Send(src, sendcnts[i], sendtype, i, MPI_ANY_TAG, comm);
         }/* end if() */
      }/* end for() */
   }/* end if() */
   else // recv data from root
   {
      MPI_Recv(recvbuf, recvcnt, recvtype, root, MPI_ANY_TAG, comm, &status);
   }/* end if() */
   return MPI_SUCCESS;
}/* end MPI_Scatterv() */

/********************************************************************
MPI_Barrier()

Synh. all ranks.
********************************************************************/
int MPI_Barrier (MPI_Comm comm)
{
   static unsigned int whichBarrier = 1;
   int i;
   int errCode = MPI_SUCCESS;
   char barrierFileName[FMPI_MAX_FNAME_SIZE];
   char barrierRoot[FMPI_MAX_FNAME_SIZE];
   FILE * pBarrier;
   unsigned int howLong;

   //master processor
   if(FMPI_gMpiRank == 0)
   {
      /*
      Wait for other processors to arrive at barrier, this is signaled
      by creation of their barrier files.
      */
      sprintf(barrierRoot, "runFileMPI.Barrier.%d.", whichBarrier);
      for(i = 1; i < FMPI_gMpiSize; i++)
      {
         sprintf(barrierFileName, "%s%d", barrierRoot, i);
         howLong = 0;
         while(howLong < FMPI_BARRIER_TIMEOUT_MS)
         {
            FMPI_CheckForAbort();

            pBarrier = fopen(barrierFileName, "r");
            if(pBarrier != NULL)
            {
               fclose(pBarrier);
               break;
            }
            else
            {
               FMPI_Sleep(FMPI_BARRIER_POLL_INTERVAL_MS);
               howLong += FMPI_BARRIER_POLL_INTERVAL_MS;
               if(howLong >= FMPI_BARRIER_TIMEOUT_MS)
               {
                  printf("Error - rank # %d timed out waiting at barrier # %d for processor # %d!\n", FMPI_gMpiRank, whichBarrier, i);
                  errCode = MPI_ERROR;
                  MPI_Abort(comm, MPI_ERROR);
               }/* end if(processor at barrier) */
            }/* end else(wait some more) */
         }/* end while(waiting at barrier) */
      }/* end for(each processor) */

      /*
      Create new barrier file, this signals to the other processors that it is
      ok to proceed past the barrier.
      */      
      sprintf(barrierFileName, "%s%d", barrierRoot, FMPI_gMpiRank);
      pBarrier = fopen(barrierFileName, "w");
      fprintf(pBarrier, "%d", FMPI_gMpiRank);
      fclose(pBarrier);      
      //printf("Created new barrier file: %s\n", barrierFileName);

      /*
      Clean up previous barrier file since we now know that no processors are waiting on it.
      */
      sprintf(barrierRoot, "runFileMPI.Barrier.%d.", whichBarrier-1);
      sprintf(barrierFileName, "%s%d", barrierRoot, FMPI_gMpiRank);
      pBarrier = fopen(barrierFileName, "r");
      if(pBarrier != NULL)
      {
         fclose(pBarrier);
         FMPI_RemoveFile(barrierFileName);
      }
      //printf("Cleaned up previous barrier file: %s\n", barrierFileName);
   }/* end master processor */
   //other processors
   else
   {
      /*
      Create a uniquely named barrier file. This lets the master processor know
      that they are waiting at the barrier.
      */
      sprintf(barrierRoot, "runFileMPI.Barrier.%d.", whichBarrier);
      sprintf(barrierFileName, "%s%d", barrierRoot, FMPI_gMpiRank);
      pBarrier = fopen(barrierFileName, "w");
      fprintf(pBarrier, "%d", FMPI_gMpiRank);
      fclose(pBarrier);

      /*
      Wait for other processors to arrive at barrier, this is signaled
      by creation of the master processor's barrier file.
      */
      sprintf(barrierFileName, "%s%d", barrierRoot, 0);
      howLong = 0;
      while(howLong < FMPI_BARRIER_TIMEOUT_MS)
      {
         FMPI_CheckForAbort();

         pBarrier = fopen(barrierFileName, "r");
         if(pBarrier != NULL)
         {
            fclose(pBarrier);
            break;
         }
         else
         {
            FMPI_Sleep(FMPI_BARRIER_POLL_INTERVAL_MS);
            howLong += FMPI_BARRIER_POLL_INTERVAL_MS;
            if(howLong >= FMPI_BARRIER_TIMEOUT_MS)
            {
               printf("Error - rank # %d timed out waiting at barrier # %d for processor # 0!\n", 
                      FMPI_gMpiRank, whichBarrier);
               errCode = MPI_ERROR;
               MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
            }/* end if(processor at barrier) */
         }/* end else(wait some more) */
      }/* end while(waiting at barrier) */

      /*
      Done waiting, go ahead and clean up the barrier file
      */
      sprintf(barrierFileName, "%s%d", barrierRoot, FMPI_gMpiRank);
      FMPI_RemoveFile(barrierFileName);
   }/* end other processors */

   whichBarrier++;
   return errCode;
}/* end MPI_Barrier() */

/********************************************************************
MPI_Bcast()

Send msg to all ranks.
********************************************************************/
int MPI_Bcast(void * buf, int count, MPI_Datatype datatype, int root, 
			  MPI_Comm comm)
{
   int i;
   char BcastFileName[FMPI_MAX_FNAME_SIZE];
   FILE * pBcast;
   int * pInts = (int *)buf;
   char * pChars = (char *)buf;
   double * pDbls = (double *)buf;

   sprintf(BcastFileName, "runFileMPI.Broadcast.%d", root);
  
   //root processor prepares data file
   if(root == FMPI_gMpiRank)
   {
      pBcast = fopen(BcastFileName, "w");
      switch(datatype)
      {
         case(MPI_INTEGER) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pBcast, "%d\n", pInts[i]);
            }/* end for() */
            break;
         }
         case(MPI_CHAR) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pBcast, "%c\n", pChars[i]);
            }/* end for() */
            break;
         }
         case(MPI_DOUBLE) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pBcast, "%.32E\n", pDbls[i]);
            }/* end for() */
            break;
         }
      }/* end switch() */
      fclose(pBcast);
   }/* end if() */

   MPI_Barrier(comm);

   //non-root processors read data file
   if(FMPI_gMpiRank != root)
   {
      pBcast = fopen(BcastFileName, "r");
      switch(datatype)
      {
         case(MPI_INTEGER) :
         {
            for(i = 0; i < count; i++)
            {
               fscanf(pBcast, "%d\n", &(pInts[i]));
            }/* end for() */
            break;
         }
         case(MPI_CHAR) :
         {
            for(i = 0; i < count; i++)
            {
               fscanf(pBcast, "%c\n", &(pChars[i]));
            }/* end for() */
            break;
         }
         case(MPI_DOUBLE) :
         {
            for(i = 0; i < count; i++)
            {
               fscanf(pBcast, "%lf\n", &(pDbls[i]));
            }/* end for() */
            break;
         }
      }/* end switch() */
      fclose(pBcast);
   }/* end if() */

   MPI_Barrier(comm);

   //root processor deletes data file
   if(FMPI_gMpiRank == root)
   {
      FMPI_RemoveFile(BcastFileName);
   }

   return MPI_SUCCESS;
}/* end MPI_Bcast() */

/********************************************************************
MPI_Allreduce()

Perform an "all" reduction operation.
********************************************************************/
int MPI_Allreduce(void * sendbuf, void * recvbuf, int count,
   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
   int i;
   for (i = 0; i < FMPI_gMpiSize; i++)
   {
      MPI_Reduce(sendbuf, recvbuf, count, datatype, op, i, comm);
   }/* end for() */
   return MPI_SUCCESS;

}/* end MPI_Allreduce() */

/********************************************************************
MPI_Reduce()

Perform a reduction operation.
********************************************************************/
int MPI_Reduce(void * sendbuf, void * recvbuf, int count, 
			   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
   int i, id, tmpInt;
   char tmpChar;
   double tmpDbl;
   char ReduceFileName[FMPI_MAX_FNAME_SIZE];
   FILE * pReduce;
   int * pSendInts = (int *)sendbuf;
   char * pSendChars = (char *)sendbuf;
   double * pSendDbls = (double *)sendbuf;
   int * pRecvInts = (int *)recvbuf;
   char * pRecvChars = (char *)recvbuf;
   double * pRecvDbls = (double *)recvbuf;

   //printf("Rank # %d is performing a reduction\n", gMpiRank);

   sprintf(ReduceFileName, "runFileMPI.Reduce.%d", FMPI_gMpiRank);

   //non-root processors prepare their portions of the data file
   if(root != FMPI_gMpiRank)
   {
      pReduce = fopen(ReduceFileName, "w");
      
      switch(datatype)
      {
         case(MPI_INTEGER) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pReduce, "%d\n", pSendInts[i]);
            }/* end for() */
            break;
         }
         case(MPI_CHAR) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pReduce, "%c\n", pSendChars[i]);
            }/* end for() */
            break;
         }
         case(MPI_DOUBLE) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pReduce, "%.16E\n", pSendDbls[i]);
            }/* end for() */
            break;
         }
      }/* end switch() */
      fclose(pReduce);
   }/* end if() */
   else /* initialize the reduction */
   {
      switch(datatype)
      {
         case(MPI_INTEGER) :
         {
            for(i = 0; i < count; i++)
            {
               pRecvInts[i] = pSendInts[i];
            }/* end for() */
            break;
         }
         case(MPI_CHAR) :
         {
            for(i = 0; i < count; i++)
            {
               pRecvChars[i] = pSendChars[i];
            }/* end for() */
            break;
         }
         case(MPI_DOUBLE) :
         {
            for(i = 0; i < count; i++)
            {
               pRecvDbls[i] = pSendDbls[i];
            }/* end for() */
            break;
         }
      }/* end switch() */
   }/* end else() */

   MPI_Barrier(comm);

   //root processor reads data files and computes the reduction
   if(FMPI_gMpiRank == root)
   {
      for(id = 0; id < FMPI_gMpiSize; id++)
      {
         if(id != FMPI_gMpiRank)
         {
            sprintf(ReduceFileName, "runFileMPI.Reduce.%d", id);
            pReduce = fopen(ReduceFileName, "r");

            switch(datatype)
            {
               case(MPI_INTEGER) :
               {
                  for(i = 0; i < count; i++)
                  {
                     fscanf(pReduce, "%d\n", &tmpInt);
                     pRecvInts[i] = FMPI_IntegerReduction(pRecvInts[i], tmpInt, op);
                  }/* end for() */
                  break;
               }
               case(MPI_CHAR) :
               {
                  for(i = 0; i < count; i++)
                  {
                     fscanf(pReduce, "%c\n", &tmpChar);
                     pRecvChars[i] = FMPI_CharacterReduction(pRecvChars[i], tmpChar, op);
                  }/* end for() */
                  break;
               }
               case(MPI_DOUBLE) :
               {
                  for(i = 0; i < count; i++)
                  {
                     fscanf(pReduce, "%lf\n", &tmpDbl);
                     pRecvDbls[i] = FMPI_DoubleReduction(pRecvDbls[i], tmpDbl, op);
                  }/* end for() */
                  break;
               }
            }/* end switch(datatype) */
            fclose(pReduce);
            FMPI_RemoveFile(ReduceFileName);
         }/* end if(not root) */
      }/* end for(each processor) */
   }/* end if(root) */

   MPI_Barrier(comm);

   return MPI_SUCCESS;
}/* end MPI_Reduce() */

/********************************************************************
MPI_Recv()

Receive a msg.
********************************************************************/
int MPI_Recv(void * buf, int count, MPI_Datatype datatype, int source, 
             int tag, MPI_Comm comm, MPI_Status * status)
{
   unsigned int msgid;
   int * sortedSources;
   int i, nTags;
   int * tagList = NULL;
   char RecvFileName[FMPI_MAX_FNAME_SIZE];
   FILE * pRecvFile = NULL;
   int errCode = MPI_SUCCESS;
   int * pRecvInts = (int *)buf;
   char * pRecvChars = (char *)buf;
   double * pRecvDbls = (double *)buf;
   int bFoundAnySourceAnyTag = 0;

   unsigned int howLong = 0;
   int srcIdIndex;
   int srcId = source;
   int tagId = tag;

   //sending processor prepares the data file
   if(source != FMPI_gMpiRank)
   {
      //printf("Rank # %d is receiving from rank # %d\n", FMPI_gMpiRank, source);
      while(howLong < FMPI_BARRIER_TIMEOUT_MS)
      {
         FMPI_CheckForAbort();

         strcpy(RecvFileName, "");

         if((source == MPI_ANY_SOURCE) && (tag == MPI_ANY_TAG))
         {
            bFoundAnySourceAnyTag = 0;

            /* ---------------------------------------------------------------------------------
            Arrange source ids in ascending order so that oldest messages are processed first.
            --------------------------------------------------------------------------------- */
            sortedSources = FMPI_SortSources(FMPI_gMsgIds, FMPI_gMpiSize, FMPI_gMpiRank);

            for(srcIdIndex = 0; srcIdIndex < FMPI_gMpiSize; srcIdIndex++)
            {
               srcId = sortedSources[srcIdIndex];

               msgid = FMPI_gMsgIds[srcId][FMPI_gMpiRank];
               nTags = FMPI_CountSendTags();
               if(nTags > 0)
               {     
                  tagList = (int *)(malloc(nTags*sizeof(int)));
                  nTags = FMPI_GetSendTagList(tagList, nTags);
                  for(tagId = 0; tagId < nTags; tagId++)
                  {
                     sprintf(RecvFileName, "runFileMPI.Send.%d.%d.%d.%d.%d", srcId, FMPI_gMpiRank, tagList[tagId], comm, msgid);
                 
                     if(FMPI_FileExists(RecvFileName) == 1)
                     {
                        bFoundAnySourceAnyTag = 1;
                        break;
                     }
                  }/* end for(each tag) */
                  free(tagList);
               }/* end if(nTags > 0) */
               if(bFoundAnySourceAnyTag == 1)
               {
                  break;
               }
            }/* end for(each source) */

            free(sortedSources);
         }/* end if(ANY_SOURCE, ANY_TAG) */
         else if((source == MPI_ANY_SOURCE) && (tag != MPI_ANY_TAG))
         {
           /* ---------------------------------------------------------------------------------
            Arrange source ids in ascending order so that oldest messages are processed first.
            --------------------------------------------------------------------------------- */
            sortedSources = FMPI_SortSources(FMPI_gMsgIds, FMPI_gMpiSize, FMPI_gMpiRank);

            for(srcIdIndex = 0; srcIdIndex < FMPI_gMpiSize; srcIdIndex++)
            {
               srcId = sortedSources[srcIdIndex];
               msgid = FMPI_gMsgIds[srcId][FMPI_gMpiRank];
               sprintf(RecvFileName, "runFileMPI.Send.%d.%d.%d.%d.%d", srcId, FMPI_gMpiRank, tag, comm, msgid);
               
               if(FMPI_FileExists(RecvFileName) == 1)
               {
                  break;
               }
            }/* end for(each source) */

            free(sortedSources);
         }/* end if(ANY_SOURCE, SPECIFIC_TAG) */
         else if((source != MPI_ANY_SOURCE) && (tag == MPI_ANY_TAG))
         {
            msgid = FMPI_gMsgIds[source][FMPI_gMpiRank];
            nTags = FMPI_CountSendTags();           
            tagList = (int *)(malloc(nTags*sizeof(int)));
            if(nTags > 0)
            {     
               nTags = FMPI_GetSendTagList(tagList, nTags);
               for(tagId = 0; tagId < nTags; tagId++)
               {
                  sprintf(RecvFileName, "runFileMPI.Send.%d.%d.%d.%d.%d", source, FMPI_gMpiRank, tagList[tagId], comm, msgid);
                  //printf("Rank %d : looking  for %s\n", gMpiRank, RecvFileName); 
               
                  if(FMPI_FileExists(RecvFileName) == 1)
                  {
                     break;
                  }
               }/* end for(each tag) */
               free(tagList);
            }/* end if(nTags > 0) */
         }/* end if(SPECIFIC_SOURCE, ANY_TAG) */
         else if((source != MPI_ANY_SOURCE) && (tag != MPI_ANY_TAG))
         {
            msgid = FMPI_gMsgIds[source][FMPI_gMpiRank];
            sprintf(RecvFileName, "runFileMPI.Send.%d.%d.%d.%d.%d", source, FMPI_gMpiRank, tag, comm, msgid);
            //printf("Rank %d : receving file %s\n", FMPI_gMpiRank, RecvFileName);
         }/* end if(SPECIFIC_SOURCE, SPECIFIC_TAG) */
         else //unknown combination!
         {
            printf("MPI_Recv() : Unknown source/tag combo\n");
            printf("source = %d\n", source);
            printf("dest   = %d\n", FMPI_gMpiRank);
            printf("tag    = %d\n", tag);
            printf("comm   = %d\n", comm);
            MPI_Abort(comm, MPI_ERROR);
         }/* end if(bad source/tag combo) */

         if(FMPI_FileExists(RecvFileName) == 1) //got the desired message file
         {
            pRecvFile = FMPI_OpenFileForRecv(RecvFileName);

            status->MPI_SOURCE = srcId;
            status->MPI_TAG = tagId;
            FMPI_gMsgIds[srcId][FMPI_gMpiRank] += 1;

            switch(datatype)
            {
               case(MPI_INTEGER) :
               {
                  for(i = 0; i < count; i++)
                  {
                     if(feof(pRecvFile))
                     {
                        printf("Rank %d : unexpected end of file (%s) in MPI_Recv()\n", FMPI_gMpiRank, RecvFileName);
                        printf("Rank %d : expected %d integers but only got %d\n", FMPI_gMpiRank, count, i);
                        fclose(pRecvFile);
                        MPI_Abort(comm, MPI_ERROR);
                     }
                     fscanf(pRecvFile, "%d\n", &(pRecvInts[i]));
                  }/* end for() */
                  break;
               }
               case(MPI_CHAR) :
               {
                  for(i = 0; i < count; i++)
                  {
                     if(feof(pRecvFile))
                     {
                        printf("Rank %d : unexpected end of file (%s) in MPI_Recv()\n", FMPI_gMpiRank, RecvFileName);
                        printf("Rank %d : expected %d characters but only got %d\n", FMPI_gMpiRank, count, i);
                        fclose(pRecvFile);
                        MPI_Abort(comm, MPI_ERROR);
                     }
                     fscanf(pRecvFile, "%c\n", &(pRecvChars[i]));
                  }/* end for() */
                  pRecvChars[i] = (char)(NULL);
                  break;
               }
               case(MPI_DOUBLE) :
               {
                  for(i = 0; i < count; i++)
                  {
                     if(feof(pRecvFile))
                     {
                        printf("Rank %d : unexpected end of file (%s) in MPI_Recv()\n", FMPI_gMpiRank, RecvFileName);
                        printf("Rank %d : expected %d doubles but only got %d\n", FMPI_gMpiRank, count, i);
                        fclose(pRecvFile);
                        MPI_Abort(comm, MPI_ERROR);
                     }
                     fscanf(pRecvFile, "%lf\n", &(pRecvDbls[i]));
                  }/* end for() */
                  break;
               }
            }/* end switch() */
            fclose(pRecvFile);

            //signal receipt by deleting the message file
            FMPI_RemoveFile(RecvFileName);

            break;
         }/* end if(received message) */
         else //message has not arrived, wait a bit
         {
            FMPI_Sleep(FMPI_BARRIER_POLL_INTERVAL_MS);
            howLong += FMPI_BARRIER_POLL_INTERVAL_MS;
            if(howLong >= FMPI_BARRIER_TIMEOUT_MS)
            {
               printf("MPI_Recv() : Rank # %d timed out waiting for message!\n", FMPI_gMpiRank);
               printf("source = %d\n", source); 
               printf("dest   = %d\n", FMPI_gMpiRank);
               printf("tag    = %d\n", tag); 
               printf("comm   = %d\n", comm);
               printf("type   = %d\n", datatype);
               if(source != MPI_ANY_SOURCE)
               {
                  printf("msg id = %d\n", FMPI_gMsgIds[source][FMPI_gMpiRank]);
               }
               errCode = MPI_ERROR;
               MPI_Abort(comm, MPI_ERROR);
            }/* end if(processor at barrier) */
         }/* end else(wait some more) */
      }/* end while(waiting for message) */
   }/* end if(receiving processor) */

   return MPI_SUCCESS;
}/* end MPI_Recv() */

/********************************************************************
MPI_Send()

Send a msg.
********************************************************************/
int MPI_Send(void * buf, int count, MPI_Datatype datatype, int dest, int tag, 
             MPI_Comm comm)
{
   unsigned int msgid;
   int i;
   char SendFileName[FMPI_MAX_FNAME_SIZE];
   FILE * pSendFile;
   int errCode = MPI_SUCCESS;
   int * pSendInts = (int *)buf;
   char * pSendChars = (char *)buf;
   double * pSendDbls = (double *)buf;
   unsigned int howLong = 0;
   FILE * pBarrier;

   //sending processor prepares the data file
   if(dest != FMPI_gMpiRank)
   {
      //printf("Rank # %d is sending to rank # %d\n", gMpiRank, dest);
      msgid = FMPI_gMsgIds[FMPI_gMpiRank][dest];
      sprintf(SendFileName, "runFileMPI.Send.%d.%d.%d.%d.%d", FMPI_gMpiRank, dest, tag, comm, msgid);
      FMPI_gMsgIds[FMPI_gMpiRank][dest] += 1; //increment msg id

      pSendFile = FMPI_OpenFileForSend(SendFileName);
      
      switch(datatype)
      {
         case(MPI_INTEGER) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pSendFile, "%d\n", pSendInts[i]);
            }/* end for() */
            break;
         }
         case(MPI_CHAR) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pSendFile, "%c\n", pSendChars[i]);
            }/* end for() */
            break;
         }
         case(MPI_DOUBLE) :
         {
            for(i = 0; i < count; i++)
            {
               fprintf(pSendFile, "%.16E\n", pSendDbls[i]);
            }/* end for() */
            break;
         }
      }/* end switch() */
      FMPI_CloseFileForSend(pSendFile, SendFileName);
   }/* end if() */

   while(howLong < FMPI_BARRIER_TIMEOUT_MS)
   {
      FMPI_CheckForAbort();

      pBarrier = fopen(SendFileName, "r");
      if(pBarrier == NULL) //receiving processor will delete once data has been extracted
      {
         break;
      }
      else
      {
         fclose(pBarrier);
         FMPI_Sleep(FMPI_BARRIER_POLL_INTERVAL_MS);
         howLong += FMPI_BARRIER_POLL_INTERVAL_MS;
         if(howLong >= FMPI_BARRIER_TIMEOUT_MS)
         {
            printf("Error - rank # %d timed out waiting at send barrier!\n", FMPI_gMpiRank);
            errCode = MPI_ERROR;
            MPI_Abort(comm, MPI_ERROR);
         }/* end if(processor at barrier) */
      }/* end else(wait some more) */
   }/* end while(waiting at barrier) */

   return MPI_SUCCESS;
}/* end MPI_Send() */

/********************************************************************
MPI_Finalize()

Shut down MPI.
********************************************************************/
int MPI_Finalize(void)
{
   int i;
   int result;
   char cmd[1000];
   char finalizeFileName[1000];
   FILE * pFile;
   unsigned int howLong;

   //printf("Rank # %d has entered MPI_Finalize() --- waiting for other processors\n", gMpiRank);
   MPI_Barrier(MPI_COMM_WORLD); //synchronize processors

   //create finalize message, master waits until all processors reach this point.
   sprintf(finalizeFileName, "FileMPI.Finalize.%d", FMPI_gMpiRank);
   pFile = fopen(finalizeFileName, "w");
   fprintf(pFile, "%d", FMPI_gMpiRank);
   fclose(pFile);

   if(FMPI_gMpiRank == 0)
   {
      remove(finalizeFileName);
      for(i = 1; i < FMPI_gMpiSize; i++)
      {
         sprintf(finalizeFileName, "FileMPI.Finalize.%d", i);
         howLong = 0;
         while(howLong < FMPI_FINALIZE_TIMEOUT_MS)
         {
            pFile = fopen(finalizeFileName, "r");
            if(pFile != NULL)
            {
               fclose(pFile);
               while((result = remove(finalizeFileName)) != 0)
               {
                  FMPI_Sleep(FMPI_FINALIZE_POLL_INTERVAL_MS);
                  howLong += FMPI_FINALIZE_POLL_INTERVAL_MS;
                  if(howLong >= FMPI_FINALIZE_TIMEOUT_MS)
                  {
                     printf("Error - rank # %d timed out on finalize waiting for file removal (%s)!\n", FMPI_gMpiRank, finalizeFileName);
                     break;
                  }/* end if(processor at barrier) */
               }
               break;
            }
            else
            {
               FMPI_Sleep(FMPI_FINALIZE_POLL_INTERVAL_MS);
               howLong += FMPI_FINALIZE_POLL_INTERVAL_MS;
               if(howLong >= FMPI_FINALIZE_TIMEOUT_MS)
               {
                  printf("Error - rank # %d timed out on finalize waiting for processor # %d!\n", FMPI_gMpiRank, i);
               }/* end if(processor at barrier) */
            }/* end else(wait some more) */
         }/* end while(waiting at barrier) */
      }/* end for(each processor) */     
      //printf("Ready to finalize!\n");
   }/* end if() */
   
   //printf("Rank # %d has exited barrier\n", gMpiRank);
   //each processor cleans up comm files
   #ifdef WIN32
      sprintf(cmd, "del runFileMPI.*.%d 1>NUL 2>NUL", FMPI_gMpiRank);
      system(cmd);
   #else
      sprintf(cmd, "rm -f runFileMPI.*.%d 1>/dev/null 2>/dev/null", FMPI_gMpiRank);
      system(cmd);
   #endif
   //printf("Rank # %d has shut down MPI communication\n", gMpiRank);
   for(i = 0; i < FMPI_gMpiSize; i++)
   {
      free(FMPI_gMsgIds[i]);
   }
   free(FMPI_gMsgIds);
   FMPI_gMpiRank = 0;
   FMPI_gMpiSize = 0;
   FMPI_gMyPid = -1;
   FMPI_gMpiIsInitialized = 0;

   remove("FileMpiRunWin32.bat");
   remove("FileMpiIn.txt");
	return 0;
}/* end MPI_Finalize() */

/********************************************************************
FMPI_CountFiles()

Count the number of files that have the matching prefix in their name.
********************************************************************/
int FMPI_CountFiles(char * prefix)
{
   int pid = FMPI_GetPid();
   int count = 0;
   char tmpFileName[FMPI_MAX_FNAME_SIZE];
   char cmdLine[FMPI_MAX_LINE_SIZE];
   char tmpStr[FMPI_MAX_LINE_SIZE];
   FILE * pOut;
   sprintf(tmpFileName, "runFileMPI.dirlist.%d", pid);
   #ifdef WIN32
     sprintf(cmdLine, "dir /OD /B %s* 1>%s 2>NUL", prefix, tmpFileName);
   #else
     sprintf(cmdLine, "ls -1tr %s* 1>%s 2>/dev/null", prefix, tmpFileName);
   #endif
   system(cmdLine);
   pOut = fopen(tmpFileName, "r");
   if(pOut == NULL) return 0;
   while(1)
   {
      FMPI_CheckForAbort();

      fgets(tmpStr, FMPI_MAX_LINE_SIZE, pOut);
      if(feof(pOut))
      {
         break;
      }
     count++;
   }
   fclose(pOut);
   FMPI_RemoveFile(tmpFileName);
   return count;
}/* end FMPI_CountFiles() */

/********************************************************************
FMPI_GetProcessList()

Create a list of process ids for each MPI application.
********************************************************************/
int FMPI_GetProcessList(unsigned int * pList, int n, char * prefix)
{
   int pid = FMPI_GetPid();
   int i, idx;
   char tmpFileName[FMPI_MAX_FNAME_SIZE];
   char cmdLine[FMPI_MAX_LINE_SIZE];
   char tmpStr[FMPI_MAX_LINE_SIZE];
   FILE * pOut;

   for(i = 0; i < n; i++) pList[i] = 0;

   /* create file containing list of pids */
   sprintf(tmpFileName, "runFileMPI.dirlist.%d", pid);
   #ifdef WIN32
      sprintf(cmdLine, "dir /OD /B %s* > %s", prefix, tmpFileName);
   #else
      sprintf(cmdLine, "ls -1tr %s* > %s", prefix, tmpFileName);
   #endif
   system(cmdLine);

   pOut = fopen(tmpFileName, "r");
   if(pOut == NULL) return 0;
   idx = strlen(prefix);
   for(i = 0; i < n; i++)
   {
      fgets(tmpStr, FMPI_MAX_LINE_SIZE, pOut);
      if(feof(pOut))
      {
         break;
      }
      pList[i] = (unsigned int)(atoi(&(tmpStr[idx])));
   }/* end for() */
   fclose(pOut);
   FMPI_RemoveFile(tmpFileName);
   return i;
}/* end FMPI_GetProcessList() */

/********************************************************************
FMPI_CountSendTags()

Count the number of tags used in send operations targeted for the 
processor. Duplicate tags will be counted, so actual number of tags 
could be less.
********************************************************************/
int FMPI_CountSendTags(void)
{
   int count = 0;
   char tmpFileName[FMPI_MAX_FNAME_SIZE];
   char cmdLine[FMPI_MAX_LINE_SIZE];
   char tmpStr[FMPI_MAX_LINE_SIZE];
   char prefix[FMPI_MAX_LINE_SIZE];
   FILE * pOut;

   sprintf(tmpFileName, "runFileMPI.sendTags.%d", FMPI_gMpiRank);
   sprintf(prefix, "runFileMPI.Send.*.%d.*.*", FMPI_gMpiRank);

   #ifdef WIN32
     sprintf(cmdLine, "dir /OD /B %s 1>%s 2>NUL", prefix, tmpFileName);
   #else
     sprintf(cmdLine, "ls -1tr %s 1>%s 2>/dev/null | sort -u", prefix, tmpFileName);
   #endif
   system(cmdLine);
   pOut = fopen(tmpFileName, "r");
   if(pOut == NULL) return 0;
   while(1)
   {
      FMPI_CheckForAbort();

      fgets(tmpStr, FMPI_MAX_LINE_SIZE, pOut);
      if(feof(pOut))
      {
         break;
      }
     count++;
   }
   fclose(pOut);
   FMPI_RemoveFile(tmpFileName);
   return count;
}/* end FMPI_CountSendTags() */

/********************************************************************
FMPI_GetSendTagList()

Get a list of unique tags used in send operations targeted for the 
processor. Duplicate tags will not be stored, so actual number of 
tags stored could be less than nTags.

Returns actual number of tags stored.
********************************************************************/
int FMPI_GetSendTagList(int * pList, int nTags)
{
   int i, j, src, dst, tag, com;
   char tmpFileName[FMPI_MAX_FNAME_SIZE];
   char cmdLine[FMPI_MAX_LINE_SIZE];
   char tmpStr[FMPI_MAX_LINE_SIZE];
   char prefix[FMPI_MAX_LINE_SIZE];
   char * pTok;
   FILE * pOut;

   for(i = 0; i < nTags; i++)
   {
      pList[i] = MPI_ANY_TAG;
   }

   sprintf(tmpFileName, "runFileMPI.sendTags.%d", FMPI_gMpiRank);
   sprintf(prefix, "runFileMPI.Send.*.%d.*.*", FMPI_gMpiRank);

   #ifdef WIN32
     sprintf(cmdLine, "dir /OD /B %s 1>%s 2>NUL", prefix, tmpFileName);
   #else
     sprintf(cmdLine, "ls -1tr %s 1>%s 2>/dev/null | sort -u", prefix, tmpFileName);
   #endif
   system(cmdLine);
   pOut = fopen(tmpFileName, "r");
   if(pOut == NULL) return 0;
   i = 0;
   while(i < nTags)
   {
      fgets(tmpStr, FMPI_MAX_LINE_SIZE, pOut);      
      if(feof(pOut))
      {
         break;
      }
      pTok = tmpStr;
      pTok += strlen("runFileMPI.Send.");
      sscanf(pTok, "%d.%d.%d.%d", &src, &dst, &tag, &com);
      //printf("rank %d : %s --> %d %d %d %d\n", gMpiRank, tmpStr, src, dst, tag, com);
      for(j = 0; j < i; j++) //unique tags only
      {
         if(tag == pList[j])
         {
            break;
         }
      }/* end for() */
      if(j == i) //no match was found
      {
         pList[i] = tag;
         i++;
      }
   }/* end while() */
   fclose(pOut);
   FMPI_RemoveFile(tmpFileName);
   return i;
}/* end FMPI_GetSendTagList() */

/********************************************************************
FMPI_OpenFileForSend()

Create a file and open it for writing. Also creates a lock file to
facilitate synchroization of reads by other ranks.
********************************************************************/
FILE * FMPI_OpenFileForSend(char * fname)
{
   FILE * pFile;
   FILE * pLock;
   char lName[FMPI_MAX_FNAME_SIZE]; //lock file

   //create a lock file
   sprintf(lName, "%s.lock", fname);
   pLock = fopen(lName, "w");
   fprintf(pLock, "rank %d : locked", FMPI_gMpiRank);
   fclose(pLock);

   //open file for writing
   pFile = fopen(fname, "w");
   if(pFile == NULL)
   {
      FMPI_RemoveFile(lName);
      return NULL;
   }
   else
   {
      return(pFile);
   }
}/* end FMPI_OpenFileForSend() */

/********************************************************************
FMPI_CloseFileForSend()

Close an open file and remove the associated lock file.
********************************************************************/
void FMPI_CloseFileForSend(FILE * pFile, char * fname)
{
   char lName[FMPI_MAX_FNAME_SIZE];

   fclose(pFile);

   sprintf(lName, "%s.lock", fname);
   FMPI_RemoveFile(lName);
}/* end FMPI_CloseFileForSend() */

/********************************************************************
FMPI_OpenFileForRecv()

Open a file containing a message sent by another rank. This function
will wait for the lock file associated with the msg to be cleared
before opening the file. This ensures entire msg is written before
attempting to read it.
********************************************************************/
FILE * FMPI_OpenFileForRecv(char * fname)
{
   unsigned int howLong;
   FILE * pFile;
   char lName[FMPI_MAX_FNAME_SIZE];

   //wait for existence of file
   howLong = 0;
   while(FMPI_FileExists(fname) == 0)
   {
      FMPI_CheckForAbort();

      FMPI_Sleep(FMPI_RECV_POLL_INTERVAL_MS);
      howLong += FMPI_RECV_POLL_INTERVAL_MS;
      if(howLong >= FMPI_RECV_TIMEOUT_MS)
      {
         printf("Error - rank # %d timed out waiting on file read!\n", FMPI_gMpiRank);
         MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
         return NULL;
      }/* end if(processor at barrier) */
   }/* end while() */

   //wait for lock file to be released
   howLong = 0;
   sprintf(lName, "%s.lock", fname);
   while(FMPI_FileExists(lName) == 1)
   {
      FMPI_CheckForAbort();

      FMPI_Sleep(FMPI_RECV_POLL_INTERVAL_MS);
      howLong += FMPI_RECV_POLL_INTERVAL_MS;
      if(howLong >= FMPI_RECV_TIMEOUT_MS)
      {
         printf("Error - rank # %d timed out waiting on lock file!\n", FMPI_gMpiRank);
         MPI_Abort(MPI_COMM_WORLD, MPI_ERROR);
         return NULL;
      }/* end if(processor at barrier) */
   }/* end while() */

   //safe to open file for reading
   pFile = fopen(fname, "r");
   return pFile;
}/* end FMPI_OpenFileForRecv() */

/********************************************************************
FMPI_FileExists()

Return 0 if file can't be opened for reading, 1 otherwise.
********************************************************************/
int FMPI_FileExists(char * fname)
{
   FILE * pFile = fopen(fname, "r");
   if(pFile == NULL) return 0;
   fclose(pFile);
   return 1;
}/* end FMPI_FileExists() */


#endif /* USE_FILE_MPI */


