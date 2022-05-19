#!/usr/bin/env bash

CHEMINP=chem_upper.inp
THERMINP=therm_upper.dat
FINALFILE=mechanism.cpp

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo "Compiling ${FINALFILE}..."
