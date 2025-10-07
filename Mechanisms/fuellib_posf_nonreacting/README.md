# FuelLib POSF Non-Reacting Mechanism
This 93 species mechanism is designed to be used with the POSF fuels available in [FuelLib](https://github.com/nrel/fuellib). It provides thermo and transport data for the specific reference compounds (see ref-compounds.csv) selected for each GCxGC bin.

The files included in `spray_input_files` are are generated using FuelLib's `Export4Pele.py` script: 

~~~
cd $FUELLIB_DIR/source
python Export4Pele.py --fuel_name posf10264 
python Export4Pele.py --fuel_name posf10264 --liq_prop_model mp
~~~

More information about generating these input files can be found in [FuelLib's documentation](https://nrel.github.io/FuelLib/tutorials-export4pele.html). Note that the units in these files are in MKS as is required for PeleLMeX. Similar files can be created for PeleC using the option `--units cgs` in the `Export4Pele.py` script. As of 10/2025, PeleC only supports the MP liquid property model. Below is an example of how to generate the required input file for use with PeleC:

~~~
cd $FUELLIB_DIR/source
python Export4Pele.py --fuel_name posf10264 --units cgs --liq_prop_model mp
~~~

