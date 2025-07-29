#!/usr/bin/env bash

MECH_HOME="$(pwd)"
MECH_FILE="${MECH_HOME}/reducedS26R134_0.yaml"
bash ../converter.sh -f "${MECH_FILE}"
