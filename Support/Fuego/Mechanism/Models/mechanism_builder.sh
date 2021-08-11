#!/usr/bin/env bash

function abspath() {
    if which realpath > /dev/null; then
        if [ -d "$1" ]; then
            echo "$(realpath "$1")";
        else
            echo "$(realpath .)";
        fi;
    else
        pushd . > /dev/null;
        if [ -d "$1" ]; then
            cd "$1" || exit;
            dirs -l +0;
        else
            cd "$(dirname "$1")" || exit;
            cur_dir=$(dirs -l +0);
            if [ "$cur_dir" == "/" ]; then
                echo "$cur_dir$(basename "$1")";
            else
                echo "$cur_dir/$(basename "$1")";
            fi;
        fi;
        popd > /dev/null || exit;
    fi;
}

if [ -z "${PELE_PHYSICS_HOME+xxx}" ]; then
    PELE_PHYSICS_HOME="$(abspath ../../../../..)"
    echo "Using PELE_PHYSICS_HOME: ${PELE_PHYSICS_HOME}. Set PELE_PHYSICS_HOME you want a different one."
else
    echo "Using PELE_PHYSICS_HOME: ${PELE_PHYSICS_HOME}"
fi

export FMC="${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py"
export HEADERDIR="${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header"

if [ -z "${FUEGO_PYTHON+xxx}" ]; then
    export FUEGO_PYTHON=$(which python)
    echo "Using environment python: $FUEGO_PYTHON. Set FUEGO_PYTHON if you want a different one."
else
    echo "Using FUEGO_PYTHON: $FUEGO_PYTHON"
fi
