#!/usr/bin/env bash

if [ -z "${PELE_PHYSICS_HOME+xxx}" ]; then
    PELE_PHYSICS_HOME="$(realpath ../../../../../)"
    echo "Using PELE_PHYSICS_HOME: ${PELE_PHYSICS_HOME}. Set PELE_PHYSICS_HOME you want a different one."
fi

export FMC="${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py"
export HEADERDIR="${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header"

if [ -z "${FUEGO_PYTHON+xxx}" ]; then
    export FUEGO_PYTHON=$(which python)
    echo "Using environment python: $FUEGO_PYTHON. Set FUEGO_PYTHON if you want a different one."
fi
