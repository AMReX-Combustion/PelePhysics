# Rough estimate of ignition time
./Pele3d.llvm.ex inputs/inputs.0d_firstpass
python computeIgnitionDelay.py -v -est -f inputs/inputs.0d_firstpass

# Refied estimate of ignition time
./Pele3d.llvm.ex inputs/inputs.0d_refine
python computeIgnitionDelay.py -ref -f inputs/inputs.0d_refine
