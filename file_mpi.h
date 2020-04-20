/******************************************************************************
File      : file_mpi.h
Author    : L. Shawn Matott
Copyright : 2013, L. Shawn Matott

Implements mpi functions using file-based approach so that one can compile 
MPI code in environments that don't have MPI libraries.

Version History
10-24-13    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef FILE_MPI_H
#define FILE_MPI_H

#define MPI_RESULTS_TAG (1) /* message contains results */
#define MPI_DATA_TAG    (2) /* message contains input/parameter data */
#define MPI_REQUEST_TAG (3) /* message is a request for work */
#define MPI_INDEX_TAG   (4) /* message is an index */
#define MPI_QUIT_TAG    (5) /* quit message */

#ifdef USE_FILE_MPI /* change preprocessor def if no file-based MPI is desired */

/* Keep C++ compilers from getting confused */
#ifdef __cplusplus
extern "C"
{
#endif
 
 typedef int MPI_Comm;
 typedef int MPI_Datatype;
 typedef int MPI_Op;

 typedef struct MPI_STATUS 
 {
	 int MPI_SOURCE;
	 int MPI_TAG;
 }MPI_Status;

 /* impoted from http://web.mit.edu/course/16/16.225/mpich-1.2.5.2/src/fortran/include/mpif.h.in */
 #define MPI_SUCCESS (0)
 #define MPI_ERROR (-1)
 #define MPI_SUM (102)
 #define MPI_MIN (103)
 #define MPI_MAX (104)
 #define MPI_DOUBLE (27) 
 #define MPI_INTEGER (28)
 #define MPI_INT (28)
 #define MPI_CHAR (1)
 #define MPI_ANY_SOURCE (-2)
 #define MPI_ANY_TAG (-1)
 #define MPI_COMM_WORLD (91)
 #define MPI_MAX_PROCESSOR_NAME (256) /* max chars in processor name */

 int MPI_Init(int * argc, char *** argv);
 int MPI_Abort(MPI_Comm comm, int errorcode);

 int MPI_Comm_size(MPI_Comm comm, int * size);
 int MPI_Comm_rank(MPI_Comm comm, int * rank);

 int MPI_Get_processor_name( char *name, int *resultlen );

 int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                void *recvbuf, int *recvcnts, int *displs, 
                MPI_Datatype recvtype, int root, MPI_Comm comm);

 int MPI_Allgatherv (void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcounts, int *displs, 
                     MPI_Datatype recvtype, MPI_Comm comm );

int MPI_Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               MPI_Comm comm);

int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm);

int MPI_Scatter(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm);

int MPI_Scatterv(void *sendbuf, int *sendcnts, int *displs, 
                MPI_Datatype sendtype, void *recvbuf, int recvcnt, 
                MPI_Datatype recvtype, int root, MPI_Comm comm);

int MPI_Barrier (MPI_Comm comm);

int MPI_Bcast(void * buf, int count, MPI_Datatype datatype, int root, 
	           MPI_Comm comm);

int MPI_Reduce(void * sendbuf, void * recvbuf, int count, 
	            MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

int MPI_Allreduce(void * sendbuf, void * recvbuf, int count,
   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int MPI_Recv(void * buf, int count, MPI_Datatype datatype, int source, 
	          int tag, MPI_Comm comm, MPI_Status * status);

int MPI_Send(void * buf, int count, MPI_Datatype datatype, int dest, int tag, 
	          MPI_Comm comm);

double MPI_Wtime(void);

int MPI_Finalize(void);

#ifdef __cplusplus
}
#endif

#else /* USE_FILE_MPI */
#ifdef USE_MPI_STUB
        #include "mpi_stub.h"
#else
        #include <mpi.h>
#endif /* USE_MPI_STUB */
#endif /* USE_FILE_MPI */
#endif /* FILE_MPI_H */
