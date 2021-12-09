ck2cti --input=dodecane_lu.inp --thermo=dodecane_lutherm.dat --transport=dodecane_lutran.dat --permissive
python -m cantera.cti2yaml dodecane_lu.cti

