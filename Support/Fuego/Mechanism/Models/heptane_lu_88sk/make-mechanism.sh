#/bin/bash

# Example script for generating Fuego code with Fortran
# (note the default pickler must be changed to "f" in 
# Support/Fuego/Pythia/pythia-0.4/packages/fuego/fuego/fmc/FMC.py

CHEMINP=chem_upper.inp
THERMINP=therm_upper.dat
FINALFILE=mechanism.cpp

FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP}  -name=${FINALFILE}

echo Compiling ${FINALFILE}...
