# Minimal Clean
if test -f "PPreaction.txt"; then
    rm PPreaction.txt
fi
if test -f "PPreaction_FDJ.txt"; then
    rm PPreaction_FDJ.txt
fi
if test -f "PPreaction_AJ.txt"; then
    rm PPreaction_AJ.txt
fi
if test -f "PPreaction_GMRES.txt"; then
    rm PPreaction_GMRES.txt
fi
if test -f "Pele3d.llvm.ex"; then
    rm Pele3d.llvm.ex
fi
if test -f "log"; then
    rm log
fi

# Compile
make COMP=llvm TPL
make COMP=llvm -j 2

# Execute 0D
./Pele3d.llvm.ex inputs/inputs.0d_cvode
mv PPreaction.txt PPreaction_FDJ.txt 
./Pele3d.llvm.ex inputs/inputs.0d_cvode_aJac
mv PPreaction.txt PPreaction_AJ.txt
./Pele3d.llvm.ex inputs/inputs.0d_cvode_GMRES
mv PPreaction.txt PPreaction_GMRES.txt

python compareResults.py



