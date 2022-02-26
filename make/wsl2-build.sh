#!/bin/bash

# This build script is for use with Ubuntu 16.04 on the Windows Subsystem for 
# Linux 2 (wsl2) platform. It will probably also work with any Linux 
# installation that has the gcc 5.x (or later) development tools installed.
# gcc 5.x also requires boost::filesystem as the experimental::filesystem
# namespace lacks "relative path" functionality. On Ubuntu, the necessary boost
# support can be installed using the following commands:
#    sudo apt-get install libboost-all-dev   (this may not be needed)
#    wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz
#    tar -xzvf boost_1_78_0.tar.gz
#    cd boost_1_78_0
#    ./bootstrap.sh --with-libraries=system,filesystem --prefix=/usr
#    sudo ./b2 install               (this will take a good bit of time!)
#    sudo ./b2 install link=static   (not sure is this is needed, but it's quick)

# GCC5 : Serial compile with GCC compiler stack
# GCC5_DBG : Serial compile with GCC compiler stack and debug symbols

make GCC5_DBG 2>&1 | tee my_make_GCC_DBG.log
if [ -f Ostrich ]; then
  mv Ostrich ../../bin/Ostrich_GCC_DBG
fi

make GCC5 2>&1 | tee my_make_GCC.log
if [ -f Ostrich ]; then
   mv Ostrich ../../bin/Ostrich_GCC
fi

