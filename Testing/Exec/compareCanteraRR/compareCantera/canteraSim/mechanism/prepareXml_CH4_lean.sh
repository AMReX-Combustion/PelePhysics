ck2cti --input=CH4_lean.inp --thermo=CH4_leantherm.dat --transport=CH4_leantran.dat --permissive
python -m cantera.cti2yaml CH4_lean.cti

