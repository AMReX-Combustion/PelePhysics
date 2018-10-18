#!/bin/sh
  export DEVELOPER=marc

  export DV_DIR=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia
  export BLD_ROOT=${DV_DIR}/builds
  export BLD_CONFIG=${DV_DIR}/config
  export PYTHIA_DIR=${DV_DIR}/pythia

  export EXPORT_ROOT=${DV_DIR}/pythia-0.4

  export PATH=$BLD_CONFIG/make:${PATH}:${EXPORT_ROOT}/bin
  export PYTHONPATH=${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver
  export LD_LIBRARY_PATH=${EXPORT_ROOT}/lib:${LD_LIBRARY_PATH}

# Python
  export PYTHON_VERSION=2.7
  export FUEGO_PYTHON=python${PYTHON_VERSION}

