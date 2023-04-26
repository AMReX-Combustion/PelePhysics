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

This section provides relevant notes on using these utilities across the Pele family of codes. There may be additional subtleties based on the code-specific implementation of these capabilities, in which case it will be necessary to refer to the documentation of the code of interest.
  
Premixed Flame Initialization
=============================

Pre-computed profiles from 1D freely propogating premixed flames are used to initialize a wrinkled `flamesheet <https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials_FlameSheet.html>`_ in PeleLMeX, among other problems. Right now, this capability is not used in PeleC, but similar code that accomplishes the same task using data files of the same format is applied in PeleC. The code has two parts, a data container defined in ``PMFData.{H,cpp}`` that loads and stores data from the pre-computed profile, and a function defined in ``PMF.H`` that, when provided this data structure and ther bounds of a cell of interest, returns temperature, velocity, and mole fractions from that location.

This code has two runtime parameters that may be set in the input file: ::

  pmf.datafile = pmf.dat
  pmf.do_cellAverage = 1

The first parameter specifies the path to the PMF data file. This file contains a two-line header followed by whitespace-delimited data columns in the order: position (cm), temperature (K), velocity (cm/s), density(g/cm3), species mole fractions. Sample files are provided in the relevant PeleLMeX (and PeleC) examples, and the procedure to generate these files with a provided script is described below. The second parameter specifies whether the PMF code does a finite volume-style integral over the querried cell (``pmf.do_cellAverage = 1``) or whether the code finds an interpolated value at the midpoint of the querried cell (``pmf.do_cellAverage = 0``)
  
Generating a PMF file
~~~~~~~~~~~~~~~~~~~~~

The script ``cantera_pmf_generator.py`` solves a 1D unstrained premixed flame using `Cantera <https://doi.org/10.5281/zenodo.6387882>`_ and saves it in the appropriate file format for use by the Pele codes. To use this script, first follow the :ref:`CEPTR instructions <sec_ceptr_software>` for setting up ``poetry``. Then return to the ``Utility/PMF/`` directory and run the script: ::

  poetry -C ../../Support/ceptr/ run python cantera_pmf_generator.py -m dodecane_lu -f NC12H26 -d 0.01 -o ./

This example computes a dodecane/air flame using the ``dodecane_lu`` mechanism on a 0.01 m domain, and saves the resulting file to the present directory. You may choose any mechanism included with PelePhysics, as the Cantera ``.yaml`` format is included for all mechanisms. If you choose a reduced mechanism with QSS species, the 1D solution will be obtained using the skeletal version of the mechanism without QSS species being remove, as Cantera does not support QSS by default. In this case, the script will also save a data file with QSS species removed (ending in ``-qssa-removed.dat``) that is suitable for use with the QSS mechanisms in Pele. A range of additional conditions may be specified at the command line. To see the full set of options, use the help feature: ::

  poetry -C ../../Support/ceptr/ run python cantera_pmf_generator.py --help

Note that when running Cantera, you may need to adjust the domain size to be appropriate for your conditions in order for the solver to converge. At times, you may also need to adjust the solver tolerances and other parameters that are specified within the script.

An additional script is provided to allow plotting of the PMF solutions. This script is used as follows (if no variable index is specified, temperature is plotted by default): ::

  poetry -C ../../Support/ceptr/ run python plotPMF.py <pmf-file> <variable-index>


Turbulent Inflows
=================

Placeholder. PelePhysics supports the capability of the flow solvers to have spatially and temporally varying inflow conditions based on precomputed turbulence data. 

Generating a turbulence file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A python script is used to generate a synthetic turbulence spectrum.

Plt File Management
===================

This code contains data structures used to handle data read from plt files that is utilized by the routines that allow the code to be restarted based on data from plt files.

Diagnostics
===========

Placeholder. Once the porting of diagnostics from PeleLMeX to PelePhysics/PeleC is complete, documentation can be added here.
