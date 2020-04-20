
all:
	-@ echo ""
	-@ echo "Must specify a target"
	-@ echo "   SER : Ostrich (serial, single processor)"
	-@ echo "   MPI : OstrichMPI (parallel, using MPI)"
	-@ echo "   ISO : build IsoFit"
	-@ echo "   GCC : GCC Ostrich (serial)"
	-@ echo "   GDB : GCC Ostrich (serial debug)" 
	-@ echo ""

# openmpi
OMPI_ISO:	*.o_OMPI_ISO
	g++ -o IsoFit *.o
	rm -f *.o

OMPI:	*.o_OMPI
	mpicxx -o OstrichMPI *.o
	rm -f *.o

OMPI_SER:	*.o_OMPI_SER
	g++ -o Ostrich *.o
	rm -f *.o

*.o_OMPI:	*.cpp *.c *.h
	mpicxx -c -O2 -DUSE_MPI *.cpp
	mpicxx -c -O2 -DUSE_MPI *.c

*.o_OMPI_SER:	*.cpp *.c *.h
	g++ -c -O2  -DUSE_MPI_STUB *.cpp
	g++ -c -O2  -DUSE_MPI_STUB *.c

*.o_OMPI_ISO:	*.cpp *.c *.h
	g++ -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.cpp
	g++ -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.c

GCC:	*.o_GCC
	g++ -o OstrichGCC *.o
	 rm -f *.o

GDB:	*.o_GDB
	g++ -g -O0 -o OstrichGDB *.o
	rm -f *.o

ISO:	*.o_ISO
	icpc -o IsoFit *.o
	rm -f *.o

MPI:	*.o_MPI
	mpiicpc -o OstrichMPI *.o
	rm -f *.o

SER:	*.o_SER
	icpc -o Ostrich *.o
	rm -f *.o

*.o_GCC:	*.cpp *.c *.h
	g++ -c -O2  -DUSE_MPI_STUB *.cpp
	g++ -c -O2  -DUSE_MPI_STUB *.c

*.o_GDB:        *.cpp *.c *.h
	g++ -c -g -O0  -DUSE_MPI_STUB *.cpp
	g++ -c -g -O0  -DUSE_MPI_STUB *.c

*.o_MPI:	*.cpp *.c *.h
	mpiicpc -c -O2 -DUSE_MPI *.cpp
	mpiicpc -c -O2 -DUSE_MPI *.c

*.o_SER:	*.cpp *.c *.h
	icpc -c -O2  -DUSE_MPI_STUB *.cpp
	icpc -c -O2  -DUSE_MPI_STUB *.c

*.o_ISO:	*.cpp *.c *.h
	icpc -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.cpp
	icpc -c -O2  -DISOFIT_BUILD -DUSE_MPI_STUB *.c

clean:
	rm -f *.o

