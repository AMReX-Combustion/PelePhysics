This code enable generating synthetic HIT files usable with PelePhysics
TurbInflow capabilities.

First generate the data using the python script:
./gen_hit_ic.py -k0 4 -N 128

To generate a synthetic HIT field discretized with 128 cells and most energetic eddies
at a wave number of 4.

Then compile the C++ executable (AMReX needed):
make

And the executable to generate the turbfile (adapt the input file to your needs):
./PeleTurb3d.gnu.ex input
