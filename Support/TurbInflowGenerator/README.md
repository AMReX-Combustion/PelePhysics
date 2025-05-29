This code enable generating usable with PelePhysics TurbInflow capabilities.
It has three modes:
 - Generation from synthetic homogeneous isotropic turbulence data,
 - Generation from a set of planes dumped from a PeleLMeX simulation using the
   PelePhysics DiagFramePlane Diagnostic capability
 - Generation from a periodic 3D plt file from a PeleLMeX simulation

See the PelePhysics documentation for more details on using this tool in all modes.
The example below described how to use it for synthetic turbulence data.

First generate the data using the python script:
./gen_hit_ic.py -k0 4 -N 128

To generate a synthetic HIT field discretized with 128 cells and most energetic eddies
at a wave number of 4.

Then compile the C++ executable (AMReX needed):
make

And the executable to generate the turbfile (adapt the input file to your needs):
./PeleTurb3d.gnu.ex input hit_file=hit_ic_4_128.dat input_ncell=128
