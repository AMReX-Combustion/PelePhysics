
CHEMINP=chem.inp
THERMINP=therm.dat
TRANINP=tran.dat
FINALFILE=heptane_3sp.c

CONVERT=${PELE_PHYSICS_HOME}/Support/Fuego/Convert/convert.exe
FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

CHEMLK=chem.asc
LOG=chem.log
TRANC=tran.c
CHEMC=chem.c
TRANLOG=tran.log

${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${CHEMC}
echo Compiling ${FINALFILE}...
echo " &files"  > model_files.dat
echo "   CHEMKIN_input = \"$CHEMINP\"" >> model_files.dat
echo "   THERMO_input = \"$THERMINP\"" >> model_files.dat
echo "   TRANLIB_input = \"$TRANINP\"" >> model_files.dat
echo "   CHEMKIN_linking_file = \"$CHEMLK\"" >> model_files.dat
echo "   TRANLIB_c_file = \"$TRANC\"" >> model_files.dat
echo "   log_file = \"$TRANLOG\"" >> model_files.dat
echo " /" >> model_files.dat
${CONVERT} model_files.dat 2>&1 >> $TRANLOG
cat ${CHEMC} ${TRANC} \
          ${HEADERDIR}/header.start\
          ${HEADERDIR}/header.mec   ${CHEMINP}\
          ${HEADERDIR}/header.therm ${THERMINP}\
          ${HEADERDIR}/header.trans ${TRANINP}\
          ${HEADERDIR}/header.end > ${FINALFILE}
rm -f ${CHEMC} ${CHEMLK} ${LOG} ${TRANC} ${TRANLOG} model_files.dat
