.. highlight:: rst

.. _sec:Support:

*******
Support
*******

The main purpose of PelePhysics is to contain a library of routines for shared use by the other Pele codes,
which are contained in the ``Source`` directory. In the ``Support`` directory, PelePhysics also contains
several stand-alone tools and scripts that are used to generate data and source code that support these purposes.
These tools include:

* ``ceptr``: Chemistry Evaluation for Pele Through Recasting, machine generation of mechanism code
* ``TurbFileHIT``: Isotropic turbulence initial conditions and inflows
* ``MechanismPAH``: Append PAH module to chemical mechanism data
* ``liquidProp``: Fits to NIST data for condensed phase species

ceptr
=====

CEPTR is used to generate the C++ mechanism files from Cantera-format mechanism data for all the chemical
mechanisms in the ``Mechanisms/`` directory. For a full description of CEPTR and how to use it, see the
:ref:`CEPTR subsection <sec:ceptr>` of the Chemistry section of this documentation.

.. _sec_turbfile:

TurbFileHIT
===========

This support code contains two separate pieces: a python script that generates a synthetic 3D isotropic turbulence velocity field and
an AMReX-based C++ code that converts that data into a format that can be used by the PelePhysics :ref:`TurbInflow Utility <sec_turbinflow>`
to provide turbulent inflow boundary conditions in either PeleC or PeleLMeX.

Generating an HIT File
~~~~~~~~~~~~~~~~~~~~~~

The ``gen_hit_ic.py`` python script creates a synthetic isotropic turbulence velocity field based on the methodology described
in the appendix of `Johnsen et al. (2010) J. Comp Phys. <http://dx.doi.org/10.1016/j.jcp.2009.10.028>`_. The resulting velocity
field is stored in files with names following the pattern ``hit_ic_k0_N.dat``, where ``k0`` is the most
energetic wave number and ``N`` is the size of the grid. The initial
condition for a grid of size N^3 is generated as follows:

1. velocity fluctuations generated on a 512^3 grid in wavenumber space
2. Coefficients associated to wavenumbers that cannot be represented on the desired grid are set to 0 (sharp wavenumber cutoff)
3. inverse Fourier transform of the velocity fluctuations (512^3 grid)
4. velocity fluctuations resampled on the desired grid (N^3)

This script accepts the following options: ::

  ./gen_hit_ic.py --help
  usage: gen_hit_ic.py [-h] [-k0 K0] [-N N] [-s SEED] [-p]

  Generate the velocity fluctuations for the HIT IC

  optional arguments:
    -h, --help            show this help message and exit
    -k0 K0                Wave number containing highest energy
    -N N                  Resolution
    -s SEED, --seed SEED  Random number generator seed
    -p, --plot            Save a plot of the x-velocit1y

Generating an initial condition file is as easy as: ::

  ./gen_hit_ic.py -N 16

.. Note::

   This script utilizes standard python libraries such as SciPy and Matplotlib. If you'd like to ensure that you are using appropriate
   versions of the dependencies, follow the :ref:`CEPTR instructions <sec_ceptr_software>` for setting up ``poetry`` to  manage
   dependencies, then run the script using ``poetry run -C ../ceptr/ python gen_hit_ic.py -N 16``.

Generating Inflow Files
~~~~~~~~~~~~~~~~~~~~~~~

The ``TurbFileHIT`` directory also contains C++ source code for a small program that reads the ``.dat`` files created using the python script
above and saves them as AMReX data structures that can be read by the TurbInflow utility, which uses Taylor's hypothesis to provide
temporally varying boundary data by marching through the 3rd dimension of the HIT files at a constant velocity.

Build options can be specified in the ``GNUmakefile``. Notably, if you did not recursively clone PelePhysics with the AMReX submodule, you
will need to specify a path to AMReX (``AMREX_HOME``). But in most cases, you should be able to directly build the utility simply by running: ::

  make -j

Then run the resulting executable using the provided input file: ::

  ./PeleTurb3d.gnu.ex input

Adjust the ``hit_file`` and ``input_ncell`` parameters in the input file to match the file name and number of cells from the previous step.
The ``urms0`` parameter is a scale factor to rescale the velocity fluctuations. ``TurbFile`` specified the name of the directory where the
output files will be saved. After successful execution, the output directory should contain two files: ``HDR`` and ``DAT``.


MechanismPAH
============

This functionality takes an existing set of mechanism, transport, and thermodynamic data files in Chemkin format and attempts to add a PAH module to it. It checks if any species are duplicated by comparing the atomic composition and the enthalpy curves. It check if any reactions are duplicated. If any species or reactions are duplicated, these are skipped. Once the new yaml file is created, ignition delay time at 1 and 20 atm is compared between the original and new mechanism to see what impact the PAH module has on the mechanism. Plot files of these values are created for you to decide if the differences are significant.

Usage
~~~~~

You need to have the NumPy, Cantera, and MatPlotLib python modules. In order to run, use the following ::

     python addPAHmech.py --mech origmech.inp --transport origtrans.dat --thermo origthermo.dat --fuelname NC10H22

where ``origmech.inp``, ``origthermo.dat``, and ``origtrans.dat`` are the initial mechanism, thermodynamic,
and transport files to have the PAH module amended to.

Disclaimer
~~~~~~~~~~

The resulting IDT should be studied to determine if the new mechanism is now compromised with the addition of the PAH module. This is left up to the user's discretion. This has only been tested on a few mechanisms and might have bugs.

liquidProp
==========

This is a python script that reads in an NIST property file for a condensed or saturated phase of a species.
Files for rho, mu, and lambda (thermal conductivity) should be provided in a single directory and
named ``rho.dat``, ``mu.dat``, and ``lambda.dat``. The usage for this script is::

  $ python propcoeff.py -h

  usage: propcoeff.py [-h] --species NC10H22 [--file_loc FILE_LOC] [--units UNITS] [--vars VARS [VARS ...]]

  options:
    -h, --help            show this help message and exit
    --species NC10H22     Species name
    --file_loc FILE_LOC   Location of data files. Files should be called rho.dat, mu.dat, and/or lambda.dat
    --units UNITS         Units, either MKS or CGS
    --vars VARS [VARS ...]
                          Which variables to fit, ex. mu lambda rho
