# Ignition delay calculation

Compute ignition delay using a 0D reactor. Ignition delay is the time to reach initial temperature + 400K


Computation of ignition delay occurs in 2 pass

1. An estimate of the ignition delay is computed

2. The estimate is refined

## Compute ignition delay 

1. Link your mechanism in `GNUmakefile` (here using `dodecane_lu`)

2. Make sure equivalence ratio calculation is correct for you mechanism in `GPU_misc.H`

3. Adjust your parameters in `inputs/inputs.0d_firstpass`

4. Execute `bash hrTest.sh` 
