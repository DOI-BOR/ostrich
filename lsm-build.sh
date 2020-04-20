#!/bin/bash

module load intel
# source /util/academic/intel/composer_xe_2013/bin/compilervars.sh intel64
module load intel-mpi

echo "Building serial application"
make clean
make SER 2>&1 | tee make-serial.log
cp Ostrich /projects/ccrstaff/lsmatott/research/bin/OstrichSerial

echo "Building parallel application"
make clean
make MPI 2>&1 | tee make-mpi.log
cp OstrichMPI /projects/ccrstaff/lsmatott/research/bin/OstrichMPI

echo "Building IsoFit application"
make clean
make ISO 2>&1 | tee make-isofit.log
cp IsoFit /projects/ccrstaff/lsmatott/research/bin/IsoFit

echo "Building GNU serial application"
make clean
make GCC 2>&1 | tee make-serial.log
cp OstrichGCC /projects/ccrstaff/lsmatott/research/bin/OstrichGCC

