#!/bin/bash
set -e

export INSTALL_PREFIX=$(pwd)/sundials/instdir
echo $INSTALL_PREFIX

export SUITESPARSE_PREFIX=$(pwd)/SuiteSparse

if [ ! -d sundials ]; then
  git clone https://github.com/LLNL/sundials
  cd sundials/
  git checkout develop
  git apply ../sundials-5.0.0-dev.1+cuda_nvector_allocators.patch
  git apply ../sundials-5.0.0-dev.1+cuda_nvector_allocators_part2.patch
else
  cd sundials/
fi


if [ ! -d $INSTALL_PREFIX ]; then
  mkdir instdir builddir
fi


cd builddir/
cmake    \
   -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}    \
   -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON    \
   -DCMAKE_C_COMPILER=$(which gcc)    \
   -DCMAKE_CXX_COMPILER=$(which g++)    \
   -DCUDA_HOST_COMPILER=$(which g++)   \
   -DEXAMPLES_INSTALL_PATH=${INSTALL_PREFIX}/examples    \
   -DCMAKE_BUILD_TYPE=Release    \
   -DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG"    \
   -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG"    \
   -DCUDA_ENABLE=ON   \
   -DMPI_ENABLE=OFF    \
   -DOPENMP_ENABLE=ON \
   -DKLU_ENABLE=ON \
   -DKLU_INCLUDE_DIR=${SUITESPARSE_PREFIX}/include \
   -DKLU_LIBRARY_DIR=${SUITESPARSE_PREFIX}/lib \
   -DEXAMPLES_ENABLE=OFF ../

make
make install

cd ../../
