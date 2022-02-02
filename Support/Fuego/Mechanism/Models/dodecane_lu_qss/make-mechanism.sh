#!/usr/bin/env bash

# Prepare the QSS mechanism
mkdir qssaTools/mechanismFiles
cp skeletal.inp qssaTools/mechanismFiles
cp therm.dat qssaTools/mechanismFiles
cp tran.dat qssaTools/mechanismFiles
cp non_qssa_list.txt qssaTools/mechanismFiles

cd qssaTools
bash transform_spec.sh
#bash transform_reac.sh
cp output/qssa_final.inp ..
cp output/therm_nostar.dat ..
cp output/tran_nostar.dat ..
cd ..



CHEMINP=qssa_final.inp
THERMINP=therm_nostar.dat
FINALFILE=mechanism.cpp

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo "Compiling ${FINALFILE}..."
