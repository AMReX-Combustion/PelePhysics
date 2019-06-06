
CHEMINP=JL4.inp
THERMINP=JL4.therm
FINALFILE=JL4.cpp

FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}
echo Compiling ${FINALFILE}...
