.. highlight:: rst

.. _sec:Utility:

*******
Utility
*******

In addition to routines for evaluating chemical reactions, transport properties, and equation of state functions, PelePhysics includes other shared utilities that are utilized by both PeleC and PeleLM(eX). These utilities include support for:

* Premixed Flame (``PMF``) initialization from precomputed 1D flame profiles
* Turbulent inflows (``TurbInflow``) on domain boundaries from saved turbulence data
* Plt file management (``PltFileManager``)
* Output of runtime ``Diagnostics``
* Basic ``Utilities``, including unit conversions

This section provides relevant notes on using these utilities across the Pele family of codes. There may be additional subtleties based on the code-specific implementation of these capabilities, in which case it will be necessary to refer to the documentation of the code of interest.

Premixed Flame Initialization
=============================

Pre-computed profiles from 1D freely propagating premixed flames are used to initialize a wrinkled `flamesheet <https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials_FlameSheet.html>`_ in PeleLMeX, among other problems. Right now, this capability is not used in PeleC, but similar code that accomplishes the same task using data files of the same format is applied in PeleC. The code has two parts, a data container defined in ``PMFData.{H,cpp}`` that loads and stores data from the pre-computed profile, and a function defined in ``PMF.H`` that, when provided this data structure and the bounds of a cell of interest, returns temperature, velocity, and mole fractions from that location.

This code has two runtime parameters that may be set in the input file: ::

  pmf.datafile = pmf.dat
  pmf.do_cellAverage = 1

The first parameter specifies the path to the PMF data file. This file contains a two-line header followed by whitespace-delimited data columns in the order: position (cm), temperature (K), velocity (cm/s), density(g/cm3), species mole fractions. Sample files are provided in the relevant PeleLMeX (and PeleC) examples, and the procedure to generate these files with a provided script is described below. The second parameter specifies whether the PMF code does a finite volume-style integral over the queried cell (``pmf.do_cellAverage = 1``) or whether the code finds an interpolated value at the midpoint of the queried cell (``pmf.do_cellAverage = 0``)

Generating a PMF file
~~~~~~~~~~~~~~~~~~~~~

The script ``cantera_pmf_generator.py`` solves a 1D unstrained premixed flame using `Cantera <https://doi.org/10.5281/zenodo.6387882>`_ and saves it in the appropriate file format for use by the Pele codes. To use this script, first follow the :ref:`CEPTR instructions <sec_ceptr_software>` for setting up ``poetry``. Then return to the ``Utility/PMF/`` directory and run the script: ::

  poetry -C ../../Support/ceptr/ run python cantera_pmf_generator.py -m dodecane_lu -f NC12H26 -d 0.01 -o ./

This example computes a dodecane/air flame using the ``dodecane_lu`` mechanism on a 0.01 m domain, and saves the resulting file to the present directory. You may choose any mechanism included with PelePhysics, as the Cantera ``.yaml`` format is included for all mechanisms. If you choose a reduced mechanism with QSS species, the 1D solution will be obtained using the skeletal version of the mechanism without QSS species being remove, as Cantera does not support QSS by default. In this case, the script will also save a data file with QSS species removed (ending in ``-qssa-removed.dat``) that is suitable for use with the QSS mechanisms in Pele. A range of additional conditions may be specified at the command line. To see the full set of options, use the help feature: ::

  poetry -C ../../Support/ceptr/ run python cantera_pmf_generator.py --help

Note that when running Cantera, you may need to adjust the domain size to be appropriate for your conditions in order for the solver to converge. At times, you may also need to adjust the solver tolerances and other parameters that are specified within the script.

An additional script is provided to allow plotting of the PMF solutions. This script is used as follows (if no variable index is specified, temperature is plotted by default): ::

  poetry -C ../../Support/ceptr/ run python plotPMF.py <pmf-file> <variable-index>

.. _sec_turbinflow:

Turbulent Inflows
=================

PelePhysics supports the capability of the flow solvers to have spatially and temporally varying inflow conditions based on precomputed turbulence data (3 components of velocity fluctuations).
A three-dimensional synthetic istropic turbulence is generated and Taylor's hypothesis is used to convert this data into two-dimensional planar data by moving through the third dimension at fixed velocity (the TurbInflow capability is currently not supported for 2D simulations). To reduce memory requirements, a fixed number of planes in the third dimension are read at a time and then replaced when exhausted. This number of planes can be specified at runtime to balance I/O and memory requirements. Multiple turbulent inflow patches can be applied on any domain boundary or on multiple domain boundaries. If multiple patches overlap on the same face, the data from each overlapping patch are superimposed.

This PelePhysics utility provides the machinery to load the data files and interpolate in space and time onto boundary patches, which may or may not cover an entire boundary. Additional code within PeleC and PeleLMeX is required to drive this functionality, and the documentation for the relevant code should be consulted to fully understand the necessary steps. Typically, the TurbInflow capability provides veloicity fluctuations, and the mean inflow velocity must be provided through another means. For each code, an example test case for the capability is provided in `Exec/RegTests/TurbInflow`. Note that the turbulence data is always a square of size 2pi and has a fluctuation velocity scale of unity. Inputs are available as part of this utility to rescale the data as needed. If differently shaped inlet patches are required, this must be done by masking undesired parts of the patch on the PeleC or PeleLMeX side of the implementation.

Generating a turbulence file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relevant data files are generated using the tools in ``Support/TurbFileHIT``
and usage instructions are available in the :ref:`documentation <sec_turbfile>` on these tools.

Input file options
~~~~~~~~~~~~~~~~~~

The input file options that drive this utility (and thus are inherited by PeleC and PeleLMeX) are
provided below. ::

  turbinflows=low high                                    # Names of injections (can provide any number)

  turbinflow.low.turb_file      = TurbFileHIT/TurbTEST    # Path to directory created in previous step
  turbinflow.low.dir            = 1                       # Boundary normal direction (0,1, or 2) for patch
  turbinflow.low.side           = "low"                   # Boundary side (low or high) for patch
  turbinflow.low.turb_scale_loc = 633.151                 # Factor by which to scale the spatial coordinate between the data file and simulation
  turbinflow.low.turb_scale_vel = 1.0                     # Factor by which to scale the velocity between the data file and simulation
  turbinflow.low.turb_center    = 0.005 0.005             # Center point where turbulence patch will be applied
  turbinflow.low.turb_conv_vel  = 5.                      # Velocity to move through the 3rd dimension to simulate time evolution
  turbinflow.low.turb_nplane    = 32                      # Number of planes to read and store at a time
  turbinflow.low.time_offset    = 0.0                     # Offset in time for reading through the 3rd dimension
  turbinflow.low.verbose        = 0                       # verbosity level

  turbinflow.high.turb_file      = TurbFileHIT/TurbTEST   # All same as above, but for second injection patch
  turbinflow.high.dir            = 1
  turbinflow.high.side           = "high"
  turbinflow.high.turb_scale_loc = 633.151
  turbinflow.high.turb_scale_vel = 1.0
  turbinflow.high.turb_center    = 0.005 0.005
  turbinflow.high.turb_conv_vel  = 5.
  turbinflow.high.turb_nplane    = 32
  turbinflow.high.time_offset    = 0.0006
  turbinflow.high.verbose        = 2

Plt File Management
===================

This code contains data structures used to handle data read from plt files that is utilized by the routines that allow the code to be restarted based on data from plt files.

Diagnostics
===========

Placeholder. Once the porting of diagnostics from PeleLMeX to PelePhysics/PeleC is complete, documentation can be added here.

Filter
======

A utility for filtering data stored in AMReX data structures. When initializing the ``Filter`` class, the filter type
and filter width to grid ratio are specified. A variety of filter types are supported:

* ``type = 0``: no filtering
* ``type = 1``: standard box filter
* ``type = 2``: standard Gaussian filter

We have also implemented a set of filters defined in Sagaut & Grohens (1999) Int. J. Num. Meth. Fluids:

* ``type = 3``: 3 point box filter approximation (Eq. 26)
* ``type = 4``: 5 point box filter approximation (Eq. 27)
* ``type = 5``: 3 point box filter optimized approximation (Table 1)
* ``type = 6``: 5 point box filter optimized approximation (Table 1)
* ``type = 7``: 3 point Gaussian filter approximation
* ``type = 8``: 5 point Gaussian filter approximation (Eq. 29)
* ``type = 9``: 3 point Gaussian filter optimized approximation (Table 1)
* ``type = 10``: 5 point Gaussian filter optimized approximation (Table 1)

.. warning:: This utility is not aware of EB or domain boundaries. If the filter stencil extends across these boundaries,
             the boundary cells are treated as if they are fluid cells.
             It is up to the user to ensure an adequate number of ghost cells in the arrays are appropriately populated,
             using the ``get_filter_ngrow()`` member function of the class to determine the required number of ghost cells.


Developing
~~~~~~~~~~

The weights for these filters are set in ``Filter.cpp``. To add a
filter type, one needs to add an enum to the ``filter_types`` and
define a corresponding ``set_NAME_weights`` function to be called at
initialization.

The application of a filter can be done on a Fab or MultiFab. The loop nesting
ordering was chosen to be performant on existing HPC architectures and
discussed in PeleC milestone reports. An example call to the filtering operation is

::

   les_filter = Filter(les_filter_type, les_filter_fgr);
   ...
   les_filter.apply_filter(bxtmp, flux[i], filtered_flux[i], Density, NUM_STATE);

The user must ensure that the correct number of grow cells is present in the Fab or MultiFab.

Utilities
=============================
The utilities namespace contains unit conversions, which are particularly useful for PeleLM(eX) users working with mixed unit systems. The following aliases, defined in a file header, enable straightforward conversions between MKS and CGS units for use in equation of state (``eos``) function calls:

::

   namespace m2c = pele::physics::utilities::mks2cgs;
   namespace c2m = pele::physics::utilities::cgs2mks;

For example, users can call ``eos.PYT2R`` as follows:

::

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real P_mean = 101325.0_rt; // Pressure in Pa (MKS)
   amrex::Real massfrac[NUM_SPECIES]; // Mass fractions tracked by PeleLM(eX)
   amrex::Real Temp = 300.0_rt; // Temp in K
   amrex::Real rho_cgs = 0.0_rt; // Density in g/cm^3 (CGS)
   
   // Calculate density in CGS units, then convert to MKS
   eos.PYT2R(m2c::P(P_mean), massfrac, Temp, rho_cgs);
   amrex::Real rho = c2m::Rho(rho_cgs); // Convert eos density to MKS