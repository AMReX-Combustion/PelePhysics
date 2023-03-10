ck2cti --input=chem-CH4-2step.inp --thermo=chem-CH4-2steptherm.dat --transport=chem-CH4-2steptran.dat --permissive
python -m cantera.cti2yaml chem-CH4-2step.cti

