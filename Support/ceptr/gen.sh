#/bin/bash

# Generate mechanisms concurrently

export PELE_PHYSICS_HOME=$(readlink -f ../..)
MODEL_DIR=${PELE_PHYSICS_HOME}/Support/Mechanism/Models

while read line; do
  echo "${line}"
  poetry run convert -f ${MODEL_DIR}/${line} &
done < ${MODEL_DIR}/list_mech
wait

while read line; do
  echo "${line}"
  poetry run qssa -fq ${MODEL_DIR}/${line} &
done < ${MODEL_DIR}/list_qss_mech
wait

while read line; do
  echo "${line}"
  poetry run convert -fq ${MODEL_DIR}/${line} &
done < ${MODEL_DIR}/list_qss_mech
wait
