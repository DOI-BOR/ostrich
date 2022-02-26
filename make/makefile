
# Ostrich files
SRC_DIR   	:= ../src
INCLUDE_DIR 	:= ../include
SRC_FILES_CPP	:= $(wildcard $(SRC_DIR)/*.cpp)

# Disk based MPI stub files
MPI_DISK_INCLUDE := ../extern/stubs/disk
MPI_DISK_SRC := ../extern/stubs/disk
DISK_FILES_C	:= $(wildcard $(MPI_DISK_SRC)/*.c)

# Memory based MPI stub files
MPI_MEM_INCLUDE := ../extern/stubs/memory
MPI_MEM_SRC := ../extern/stubs/memory
MEM_FILES_C	:= $(wildcard $(MPI_MEM_SRC)/*.c)


all:
	-@ echo ""
	-@ echo "Must specify a target"
	-@ echo "   ICC : Serial compile with Intel compiler stack"
	-@ echo "   GCC : Serial compile with GCC 8.x compiler stack"
	-@ echo "   GCC_DBG : Serial compile with GCC 8.x compiler stack and debug symbols"
	-@ echo "   GCC5 : Serial compile with GCC 5.x compiler stack"
	-@ echo "   GCC5_DBG : Serial compile with GCC 5.x compiler stack and debug symbols"
	-@ echo "   ICC_MPI : MPI compile with Intel compiler stack"
	-@ echo "   GCC_MPI : MPI compile with GCC compiler stack"
	-@ echo "   GCC_DB : Serial compile with GCC 8.x compiler stack with disk based MPI"
	-@ echo "   ICC_DB : Serial compile with Intel compiler stack with disk based MPI"
	-@ echo ""

GCC_DBG:    $(SRC_FILES_CPP) $(SRC_FILES_C)
	g++ -g -std=c++17 -I$(INCLUDE_DIR) -I$(MPI_MEM_INCLUDE) $(SRC_FILES_CPP) $(MEM_FILES_C) -ldl -lstdc++fs -o Ostrich
	rm -f *.o

GCC5_DBG:    $(SRC_FILES_CPP) $(SRC_FILES_C)
	g++ -g -std=c++17 -DGCC5X -I$(INCLUDE_DIR) -I$(MPI_MEM_INCLUDE) $(SRC_FILES_CPP) $(MEM_FILES_C) -ldl -L/usr/lib -lboost_system -lboost_filesystem -o Ostrich
	rm -f *.o

ICC_MPI:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	mpiicpc -std=c++17 -I$(INCLUDE_DIR) $(SRC_FILES_CPP) -ldl -lstdc++fs -o OstrichMPI
	rm -f *.o

GCC_MPI:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	mpigxx -std=c++17 -I$(INCLUDE_DIR) $(SRC_FILES_CPP) -ldl -lstdc++fs -o OstrichMPI 
	rm -f *.o

ICC:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	icpc  -std=c++17 -I$(INCLUDE_DIR) -I$(MPI_MEM_INCLUDE) $(SRC_FILES_CPP) $(MEM_FILES_C) -ldl -lstdc++fs -o Ostrich
	rm -f *.o

GCC:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	g++ -std=c++17 -I$(INCLUDE_DIR) -I$(MPI_MEM_INCLUDE) $(SRC_FILES_CPP) $(MEM_FILES_C) -ldl -lstdc++fs -o Ostrich
	rm -f *.o

GCC5:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	g++ -std=c++17 -DGCC5X -I$(INCLUDE_DIR) -I$(MPI_MEM_INCLUDE) $(SRC_FILES_CPP) $(MEM_FILES_C) -ldl -L/usr/lib -lboost_system -lboost_filesystem -o Ostrich
	rm -f *.o

GCC_DB:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	g++ -std=c++17 -I$(INCLUDE_DIR) -I$(MPI_DISK_INCLUDE) $(SRC_FILES_CPP) $(DISK_FILES_C) -ldl -lstdc++fs -o OstrichDB
	rm -f *.o

ICC_DB:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	icpc -std=c++17 -I$(INCLUDE_DIR) -I$(MPI_DISK_INCLUDE) $(SRC_FILES_CPP) $(DISK_FILES_C) -ldl -lstdc++fs -o OstrichDB
	rm -f *.o

clean:
	rm -f *.o

