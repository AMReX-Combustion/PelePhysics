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
./Pele3d.llvm.ex inputs/inputs.0d_firstpass


