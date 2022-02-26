#!/usr/bin/env bash

# Prepare the QSS mechanism
mkdir qssaTools/mechanismFiles
cp skeletal.inp qssaTools/mechanismFiles
cp therm.dat qssaTools/mechanismFiles
cp tran.dat qssaTools/mechanismFiles
cp non_qssa_list.txt qssaTools/mechanismFiles
#cp non_qssa_list_save.txt qssaTools/mechanismFiles/non_qssa_list.txt

cd qssaTools
bash transform_spec.sh
#bash transform_reac.sh
#bash transform_reac_final.sh
#cp output/reac_forward_to_remove ..
cp output/qssa_final.inp ..
cp output/therm_nostar.dat ..
cp output/tran_nostar.dat ..
cd ..



CHEMINP=qssa_final.inp
THERMINP=therm_nostar.dat
FINALFILE=mechanism.cpp
FORWARD_REAC_TO_REMOVE=reac_forward_to_remove

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE} ${FORWARD_REAC_TO_REMOVE}

echo "Compiling ${FINALFILE}..."
