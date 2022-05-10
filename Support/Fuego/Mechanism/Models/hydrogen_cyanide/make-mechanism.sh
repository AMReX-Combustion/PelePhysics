#!/usr/bin/env bash


ORIG_CHEMINP=orig/mechanism.inp
ORIG_THERMINP=orig/therm.dat
ORIG_TRANINP=orig/tran.dat

CHEMINP=mechanism.inp
THERMINP=therm.dat
FINALFILE=mechanism.cpp

cp ${ORIG_CHEMINP} ${CHEMINP}

python ../useful_script/script_thermo.py ${CHEMINP} ${ORIG_THERMINP}
mv THERMOFILE_APPEND.txt ${THERMINP}

python ../useful_script/script_trans.py ${CHEMINP} ${ORIG_TRANINP}
cat TRANFILE_APPEND.txt >> ${CHEMINP}
rm TRANFILE_APPEND.txt

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo "Compiling ${FINALFILE}..."
