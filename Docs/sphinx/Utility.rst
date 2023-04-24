.. highlight:: rst

.. _sec:Utility:

*******
Utility
*******

.. note:: Placeholder, to be written

In addition to routines for evaluating chemical reactions, transport properties, and equation of state functions, PelePhysics includes other shared utilities that are utilized by both PeleC and PeleLM(eX). These utilities include support for:

* Premixed Flame (``PMF``) initialization from precomputed 1D flame profiles
* Turbulent inflows (``TurbInflow``) on domain boundaries from saved turbulence data
* Plt file management (``PltFileManager``)
* Output of runtime ``Diagnostics``

This section provides relevant notes on using these utilities across the Pele family of codes. There may be additional subtleties based on the code-specific implementation of these capabilities, in which case it will be necessary to refer to the documentation of the code of interest.
  
Premixed Flame Initialization
=============================

Pre-computed profiles from 1D freely propogating premixed flames are used to initialize a wrinkled flame profile in the PMF regtest in PeleC and in the FlameSheet regtest/tutorial in PeleLMeX, among other problems.  

Turbulent Inflows
=================

Plt File Management
===================

Diagnostics
===========
