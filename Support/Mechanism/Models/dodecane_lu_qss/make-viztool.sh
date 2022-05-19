#!/usr/bin/env bash
export DV_DIR=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia
export EXPORT_ROOT=${DV_DIR}/pythia-0.5
export PYTHONPATH=${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver:${EXPORT_ROOT}/packages/fuego/fuego/fmc

# Prepare the QSS mechanism
qssaToolsDir=../useful_script/qssaTools
skeletalMechanism=skeletal.inp
non_qssa_list=non_qssa_list.txt
thermoFile=therm.dat
tranFile=tran.dat

# Make graph
# Write QSSA inp
python $qssaToolsDir/makeQSSAInp/make_inp_file_qssa.py -sk $skeletalMechanism -nqss $non_qssa_list
# Write remove star names
python $qssaToolsDir/makeQSSAInp/removeStarFromSpecies.py -sk $skeletalMechanism -th $thermoFile -tr $tranFile -nqss $non_qssa_list
# Make the directed graph
python $qssaToolsDir/vizTool/makeDiGraph.py -v 
python $qssaToolsDir/vizTool/makeQSSQuadGraph.py -v
echo "Made QSS Directed Graph"


