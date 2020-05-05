#!/bin/bash
set -e

export INSTALL_PREFIX=$(pwd)/sundials/instdir
echo $INSTALL_PREFIX

##export SUITESPARSE_PREFIX=$(pwd)/SuiteSparse

if [ ! -d sundials ]; then
  git clone https://github.com/LLNL/sundials
  cd sundials/
  git checkout develop
else
  cd sundials/
fi


if [ ! -d $INSTALL_PREFIX ]; then
  mkdir instdir builddir
fi


cd builddir/
#note that these compil flags have been tested on CORI with a gnu compiler 
cmake    \
   -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}    \
   -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON    \
   -DCMAKE_C_COMPILER=$(which gcc)    \
   -DCMAKE_CXX_COMPILER=$(which g++)    \
   -DCUDA_HOST_COMPILER=$(which g++)   \
   -DCMAKE_CUDA_COMPILER=$(which nvcc)   \
   -DSUNDIALS_INDEX_SIZE=32 \
   -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples    \
   -DCMAKE_BUILD_TYPE=Release    \
   -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG"    \
   -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"    \
   -DCUDA_ENABLE=ON   \
   -DMPI_ENABLE=OFF    \
   -DOPENMP_ENABLE=ON \
   -DKLU_ENABLE=ON \
   -DKLU_INCLUDE_DIR=${SUITESPARSE_DIR}/include \
   -DKLU_LIBRARY_DIR=${SUITESPARSE_DIR}/lib \
   -DEXAMPLES_ENABLE=ON ../

make
make install

cd ../../
