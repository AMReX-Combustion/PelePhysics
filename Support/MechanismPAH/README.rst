PAH Mechanism Modifier
----------------------

This functionality takes an existing set of mechanism, transport, and thermodynamic data files in Chemkin format and attempts to add a PAH module to it. It checks if any species are duplicated by comparing the atomic composition and the enthalpy curves. It check if any reactions are duplicated. If any species or reactions are duplicated, these are skipped. Once the new yaml file is created, ignition delay time at 1 and 20 atm is compared between the original and new mechanism to see what impact the PAH module has on the mechanism. Plot files of these values are created for you to decide if the differences are significant.

Usage
~~~~~

You need to have the NumPy, Cantera, and MatPlotLib python modules. In order to run, use the following ::

     python addPAHmech.py --mech origmech.inp --transport origtrans.dat --thermo origthermo.dat --fuelname NC10H22

where `origmech.inp`, `origthermo.dat`, and `origtrans.dat` are the initial mechanism, thermodynamic, and transport files to have the PAH module amended to.

Disclaimer
~~~~~~~~~~

The resulting IDT should be studied to determine if the new mechanism is now compromised with the addition of the PAH module. This is left up to the user's discretion. This has only been tested on a few mechanisms and might have bugs.
