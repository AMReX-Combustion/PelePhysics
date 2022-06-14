# Minimal Clean
if test -f "Pele3d.llvm.ex"; then
    rm Pele3d.llvm.ex
fi
if test -f "logBase.txt"; then
    rm logBase.txt
fi
if test -f "logSym.txt"; then
    rm logSym.txt
fi

# Compile
make COMP=llvm TPL
make COMP=llvm -j 2

# Execute 0D
./Pele3d.llvm.ex 
