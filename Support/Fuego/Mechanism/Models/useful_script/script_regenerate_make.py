import numpy as np
import sys
import csv

####################################################################
#
#  USAGE: "script_regenerate_make.py y/n"
#         It will return a make-mechanism.sh and a Make.package 
#         The y/n is here to generate a Make.package including the 
#         cpp routine for the GPU custom direct solver
#
####################################################################

def print_make_fuego(outfile):

    " print the shell script to generate a cpp mechanism with Fuego "
    " the name of this script is outfile                            "

    f = open(outfile, 'w')
    f.write("CHEMINP=mechanism.inp\n")
    f.write("THERMINP=therm.dat\n")
    f.write("FINALFILE=mechanism.cpp\n\n")
    f.write("FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py\n")
    f.write("HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header\n\n")
    f.write("${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}\n\n")
    f.write("echo Compiling ${FINALFILE}...\n")
    f.closed

def print_Make(outfile, cppsolver):

    " print the Make.package "

    f = open(outfile, 'w')
    f.write("CEXE_headers+=chemistry_file.H\n\n")
    if (cppsolver=="y"):
        print("Using solver.cpp")
        f.write("CEXE_sources+=mechanism.cpp solver.cpp solver_simplified.cpp\n\n")
    else:
        f.write("CEXE_sources+=mechanism.cpp\n\n")
    f.write("LIBRARIES +=\n\n")
    f.write("INCLUDE_LOCATIONS += $(CHEM_DIR)\n\n")
    f.write("VPATH_LOCATIONS +=  $(CHEM_DIR)/PMFs\n")

print_make_fuego("make-mechanism.sh")
print_Make("Make.package", sys.argv[1])
