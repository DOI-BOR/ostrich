#!/bin/bash

module load intel/20.2
module load intel-mpi
module load gcc/8.3.0

# ICC : Serial compile with Intel compiler stack
# GCC : Serial compile with GCC compiler stack
# ICC_MPI : MPI compile with Intel compiler stack
# GCC_MPI : MPI compile with GCC compiler stack
# GCC_DB : Serial compile with GCC compiler stack with disk based MPI
# ICC_DB : Serial compile with Intel compiler stack with disk based MPI

make GCC_DBG 2>&1 | tee my_make_GCC_DBG.log
if [ -f Ostrich ]; then
  mv Ostrich ../../bin/Ostrich_GCC_DBG
fi

make ICC 2>&1 | tee my_make_ICC.log
if [ -f Ostrich ]; then
  mv Ostrich ../../bin/Ostrich_ICC
fi

make GCC 2>&1 | tee my_make_GCC.log
if [ -f Ostrich ]; then
   mv Ostrich ../../bin/Ostrich_GCC
fi

make ICC_MPI 2>&1 | tee my_make_ICC_MPI.log
if [ -f OstrichMPI ]; then
   mv OstrichMPI ../../bin/Ostrich_ICC_IMPI
fi

make GCC_MPI 2>&1 | tee my_make_GCC_MPI.log
if [ -f OstrichMPI ]; then
   mv OstrichMPI ../../bin/Ostrich_GCC_IMPI
fi

