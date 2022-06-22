bash reenableAllProfile.sh
rm ceptr_main.py.lprof

poetry run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml
kernprof -l  ceptr_main.py -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa.yaml
