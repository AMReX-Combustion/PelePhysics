bash disableAllProfile.sh

poetry run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml

poetry run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa.yaml --hformat readable --remove_1 --remove_pow2 --min_op_count 3
