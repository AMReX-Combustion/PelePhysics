# Minimal Clean
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
./Pele3d.llvm.ex > log
