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
if compgen -G "Pele*.ex" > /dev/null; then
    rm Pele*.ex
fi

# Compile
make COMP=llvm TPL
make COMP=llvm -j 2


# Find first executable name available
execname=`find . -name "Pele*.ex" | head -1`

# Execute 0D
$execname inputs/inputs.0d_cvode
mv PPreaction.txt PPreaction_FDJ.txt 
$execname inputs/inputs.0d_cvode_aJac
mv PPreaction.txt PPreaction_AJ.txt
$execname inputs/inputs.0d_cvode_GMRES
mv PPreaction.txt PPreaction_GMRES.txt

python compareResults.py



