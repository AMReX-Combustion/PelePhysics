# Minimal Clean
if test -f "PPreaction.txt"; then
    rm PPreaction.txt
fi
if test -f "Pele3d.llvm.ex"; then
    rm Pele3d.llvm.ex
fi

# Compile
make COMP=llvm TPL
make COMP=llvm

# Execute 0D
#./Pele3d.llvm.ex inputs/inputsRK64
#./Pele3d.llvm.ex inputs/inputs.0d_denseDirect
#./Pele3d.llvm.ex inputs/inputs.0d_gmres
./Pele3d.llvm.ex inputs/inputs.0d_AJ

# Compare to Cantera
if test -f "compareCantera/canteraSim/output/CanteraReaction_dodecane.txt"; then
    cd compareCantera
    python compareResults.py
    cd ..
else
    source ~/miniconda/bin/activate cantera
    cd compareCantera/canteraSim
    python homoReact_dodecane.py  
    cd ../
    python compareResults.py
    cd ..
fi
