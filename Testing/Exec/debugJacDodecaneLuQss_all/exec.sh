# Minimal Clean
if test -f "logBase.txt"; then
    rm logBase.txt
fi
if test -f "logSym.txt"; then
    rm logSym.txt
fi
if test -f "Pele*.ex"; then
    rm Pele3d.*
fi

# Compile
make COMP=llvm TPL
make COMP=llvm -j 2


# Find first executable name available
execname=`find . -name "Pele*.ex" | head -1`

# Execute 0D
$execname
