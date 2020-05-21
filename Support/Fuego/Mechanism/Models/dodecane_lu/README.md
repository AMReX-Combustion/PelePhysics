# Analytically Reduced Mechanisms
!!This mechanism is the skeletal version. !!

There is a fortran routine to reduce it even more to an ARC. However in that case the production rate 
of the species should NOT be evaluated via the reactions provided in the chemkin reaction data file.

Note that if a .cpp file is generated at this point, 
it will implement the skeletal version of the reduced mechanism and it will not use the *.f routine.
