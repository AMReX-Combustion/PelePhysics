generateQSSOnly () {

    echo Taking care of skeletal at $1    

    #poetry run ck2yaml --input=$1/mechanism.inp --thermo=$1/therm.dat --transport=$1/tran.dat --permissive
    
    #poetry run convert -f $1/mechanism.yaml
    
    echo Taking care of QSS at $2   

    cp $1/mechanism.yaml $2/skeletal.yaml
    
    poetry run qssa -f $2/skeletal.yaml -n $2/non_qssa_list.yaml
    
    poetry run convert -f $2/qssa.yaml --qss_format_input $2/qssa_input.toml --qss_symbolic_jacobian
}
generate () {

    echo Taking care of skeletal at $1    

    poetry run ck2yaml --input=$1/mechanism.inp --thermo=$1/therm.dat --transport=$1/tran.dat --permissive
    
    poetry run convert -f $1/mechanism.yaml
    
    echo Taking care of QSS at $2   

    cp $1/mechanism.yaml $2/skeletal.yaml
    
    poetry run qssa -f $2/skeletal.yaml -n $2/non_qssa_list.yaml
    
    poetry run convert -f $2/qssa.yaml --qss_format_input $2/qssa_input.toml --qss_symbolic_jacobian
}

ModelRoot=${PELE_PHYSICS_HOME}/Support/Mechanism/Models

PATH_TO_CHEMKIN_DIR=${ModelRoot}/CH4_lean 
PATH_TO_CHEMKIN_DIR_QSS=${ModelRoot}/CH4_lean_qss
generateQSSOnly $PATH_TO_CHEMKIN_DIR $PATH_TO_CHEMKIN_DIR_QSS 

PATH_TO_CHEMKIN_DIR=${ModelRoot}/C1-C2-NO
PATH_TO_CHEMKIN_DIR_QSS=${ModelRoot}/C1-C2-NO_qss
generateQSSOnly $PATH_TO_CHEMKIN_DIR $PATH_TO_CHEMKIN_DIR_QSS

PATH_TO_CHEMKIN_DIR=${ModelRoot}/LuEthylene
PATH_TO_CHEMKIN_DIR_QSS=${ModelRoot}/LuEthylene_qss
generateQSSOnly $PATH_TO_CHEMKIN_DIR $PATH_TO_CHEMKIN_DIR_QSS

PATH_TO_CHEMKIN_DIR=${ModelRoot}/heptane_lu_88sk
PATH_TO_CHEMKIN_DIR_QSS=${ModelRoot}/heptane_lu_qss
generateQSSOnly $PATH_TO_CHEMKIN_DIR $PATH_TO_CHEMKIN_DIR_QSS

PATH_TO_CHEMKIN_DIR=${ModelRoot}/dodecane_lu
PATH_TO_CHEMKIN_DIR_QSS=${ModelRoot}/dodecane_lu_qss
generateQSSOnly $PATH_TO_CHEMKIN_DIR $PATH_TO_CHEMKIN_DIR_QSS
