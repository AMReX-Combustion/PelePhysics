This is a non-reactive mechanism for A1 (POSF10264), A2 (POSF10325), A3 (POSF10289), C1 (POSF11498), and air (O2, N2) with CO2 and H2O. The mechanism is based on the HyChem mechanism without low temperature chemistry available here: 
[https://web.stanford.edu/group/haiwanglab/HyChem/pages/download.html](https://web.stanford.edu/group/haiwanglab/HyChem/pages/download.html)

See the comments at the top of mechanism.inp for more details and citation
information.

The files included in `spray_input_files` are generated using FuelLib's `Export4Pele.py` script for use with the spray module of PelePhysics. More information about generating these input files can be found in [FuelLib's documentation](https://nrel.github.io/FuelLib/tutorials-export4pele.html). Note that the units in these files are in MKS as is required for PeleLMeX. Similar files can be created for PeleC using the option `--units cgs` in the `Export4Pele.py` script. As of 10/2025, PeleC only supports the MP liquid property model. Below is an example of how to generate the required input file for use with PeleC:

~~~
cd $FUELLIB_DIR/source
python Export4Pele.py --fuel_name posf10264 --units cgs --export_mix True --liq_prop_model mp
~~~