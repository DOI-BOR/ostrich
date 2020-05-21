/******************************************************************************
File      : mpi_stub.c
Author    : L. Shawn Matott
Copyright : 2002, L. Shawn Matott

Stubs out mpi functions, so that one can compile the Ostrich code without using
the MPI libraries.

Version History
11-18-02    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
******************************************************************************/
#include "mpi_stub.h"
#include <stdlib.h>

#ifdef USE_MPI_STUB /* change preprocessor def if no stubbing desired */

int MPI_Init(int * argc, char *** argv)
{
   return 0;
}

int MPI_Abort(MPI_Comm comm, int errorcode)
{
   return 0;
}

int MPI_Comm_size(MPI_Comm comm, int * size)
{
   *size = 1;
   return 0;
}

int MPI_Comm_rank(MPI_Comm comm, int * rank)
{
   *rank = 0;
   return 0;
}

int MPI_Allgatherv (void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcounts, int *displs, 
                     MPI_Datatype recvtype, MPI_Comm comm )
{
   return 0;
}

int MPI_Gatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                void *recvbuf, int *recvcnts, int *displs, 
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
   return 0;
}

int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
               void *recvbuf, int recvcnt, MPI_Datatype recvtype, 
               int root, MPI_Comm comm)
{
   return 0;
}

int MPI_Barrier (MPI_Comm comm)
{
   return 0;
}

int MPI_Bcast(void * buf, int count, MPI_Datatype datatype, int root, 
			  MPI_Comm comm)
{
   return 0;
}

int MPI_Reduce(void * sendbuf, void * recvbuf, int count, 
			   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
	return 0;
}

int MPI_Allreduce(void * sendbuf, void * recvbuf, int count,
   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
   return 0;
}

int MPI_Recv(void * buf, int count, MPI_Datatype datatype, int source, 
             int tag, MPI_Comm comm, MPI_Status * status)
{
	 return 0;
}

int MPI_Send(void * buf, int count, MPI_Datatype datatype, int dest, int tag, 
             MPI_Comm comm)
{
	 return 0;
}

int MPI_Finalize(void)
{
	return 0;
}

#endif /* USE_MPI_STUB */


