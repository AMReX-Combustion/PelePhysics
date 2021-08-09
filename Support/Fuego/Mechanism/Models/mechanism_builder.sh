#!/usr/bin/env bash

if [ -z "${PELE_PHYSICS_HOME+xxx}" ]; then
    PELE_PHYSICS_HOME="$(realpath ../../../../../)"
    echo "PELE_PHYSICS_HOME is not set, defaulting to ${PELE_PHYSICS_HOME}"
fi

export FMC="${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py"
export HEADERDIR="${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header"

if [ -z "${FUEGO_PYTHON+xxx}" ]; then
    export FUEGO_PYTHON=$(which python)
    echo "FUEGO_PYTHON is not set, using environment python: ${FUEGO_PYTHON}"
fi
