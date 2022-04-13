#!/usr/bin/env bash
export DV_DIR=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia
export EXPORT_ROOT=${DV_DIR}/pythia-0.5
export PYTHONPATH=${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver

# Prepare the QSS mechanism
qssaToolsDir=../useful_script/qssaTools
skeletalMechanism=skeletal.inp
non_qssa_list=non_qssa_list.txt
thermoFile=therm.dat
tranFile=tran.dat

# Write QSSA inp
python $qssaToolsDir/scripts/makeQSSAInp/make_inp_file_qssa.py -sk $skeletalMechanism -nqss $non_qssa_list
# Write remove star names
python $qssaToolsDir/scripts/makeQSSAInp/removeStarFromSpecies.py -sk $skeletalMechanism -th $thermoFile -tr $tranFile -nqss $non_qssa_list

## ~~~~ Method 1 : Remove enough species to remove quadratic coupling
#python $qssaToolsDir/scripts/removeQuadratic_species/removeQuadratic_method1.py
# ~~~~ Method 2 : Remove all reactions (forward and backward) that generate quadratic coupling
#python $qssaToolsDir/scripts/removeQuadratic_reactions/removeQuadratic_method2.py
# ~~~~ Method 3 : Remove all reactions (forward, backward or both) that generate quadratic coupling
python $qssaToolsDir/scripts/removeQuadratic_reactions/removeQuadratic_method3.py
cp output/reac_forward_to_remove .

cp output/qssa_final.inp .
cp output/therm_nostar.dat .
cp output/tran_nostar.dat .

# Compile mechanism

CHEMINP=qssa_final.inp
THERMINP=therm_nostar.dat
FINALFILE=mechanism.cpp

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo "Compiling ${FINALFILE}..."
