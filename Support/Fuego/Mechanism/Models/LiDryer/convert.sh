#!/usr/bin/env bash

MECH_HOME="$(pwd)"
export MECH_FILE="${MECH_HOME}/mechanism.yaml"
bash ../converter.sh
