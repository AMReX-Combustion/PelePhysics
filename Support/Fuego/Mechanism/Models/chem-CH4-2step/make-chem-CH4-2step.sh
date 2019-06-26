
CHEMINP=chem-CH4-2step.mec
THERMINP=chem-CH4-2step.therm
FINALFILE=chem-CH4-2step.cpp

FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}
echo Compiling ${FINALFILE}...
