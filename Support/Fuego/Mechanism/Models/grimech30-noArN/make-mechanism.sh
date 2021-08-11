#!/usr/bin/env bash

CHEMINP=mechanism.inp
THERMINP=therm.dat
FINALFILE=mechanism.cpp

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo "Compiling ${FINALFILE}..."
