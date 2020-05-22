
SRC_DIR   	:= ../src
INCLUDE_DIR 	:= ../include

SRC_FILES_CPP	:= $(wildcard $(SRC_DIR)/*.cpp)
SRC_FILES_C	:= $(wildcard $(SRC_DIR)/*.c)

all:
	-@ echo ""
	-@ echo "Must specify a target"
	-@ echo "   ICC : Intel compiler stack"
	-@ echo "   GCC : GCC compiler stack"
	-@ echo ""


ICC:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	mpiicpc -I$(INCLUDE_DIR) $(SRC_FILES_CPP) $(SRC_FILES_C) -o OstrichMPI
	rm -f *.o

GCC:	$(SRC_FILES_CPP) $(SRC_FILES_C)
	mpic++ -I$(INCLUDE_DIR) $(SRC_FILES_CPP) $(SRC_FILES_C) -o OstrichMPI
	rm -f *.o


clean:
	rm -f *.o

