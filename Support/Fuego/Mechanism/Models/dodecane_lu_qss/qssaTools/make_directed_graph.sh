#!/usr/bin/env bash
export DV_DIR=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia
export BLD_ROOT=${DV_DIR}/builds
export BLD_CONFIG=${DV_DIR}/config
export PYTHIA_DIR=${DV_DIR}/pythia
export EXPORT_ROOT=${DV_DIR}/pythia-0.5
export PATH=$BLD_CONFIG/make:${PATH}:${EXPORT_ROOT}/bin
export PYTHONPATH=${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver
export LD_LIBRARY_PATH=${EXPORT_ROOT}/lib:${LD_LIBRARY_PATH}


python -c "import networkx"
notAvail=$?

if [ $notAvail -eq 1 ]
then
    echo "networkx module unavailable, cannot plot directed graphs"
elif [ $notAvail -eq 0 ]
then
    # Write QSSA inp
    python scripts/makeQSSAInp/make_inp_file_qssa.py scripts/inputs/input
    # Write remove star names
    python scripts/makeQSSAInp/removeStarFromSpecies.py scripts/inputs/input_star
    # Make the directed graph
    python scripts/vizTool/makeDiGraph.py scripts/inputs/input_diGraph
    python scripts/vizTool/makeQSSQuadGraph.py scripts/inputs/input_diGraph
fi
