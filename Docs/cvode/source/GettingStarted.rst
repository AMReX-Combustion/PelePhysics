.. highlight:: rst

`PelePhysics` Quickstart
============================

Greetings impatient user. Once again, note that this documentation focuses on the CHEMISTRY part of `PelePhysics`. 
If you are already familiar with this code, and would like to know which chemistry-related keywords 
and/or environment variables to set in your various input files to perform a simulation with one of the codes 
available in the `PeleSuite`, then I invite you to directly skip to section :ref:`sec:subsPPOptions`. 
If you are in a hurry but still would like more context, visit section :ref:`sec:subsWD`
to be referred to portions of this document that are of interest to you.


For the sake of clarity, the following provides a quick reminder on how to obtain `PelePhysics` and all the supporting software required. A more comprehensive `PelePhysics` manual is on the way.


Obtaining `PelePhysics`
----------------------------

First, make sure that "git" is installed on your machine---we recommend version 1.7.x or higher. Then...

1. Download the `AMReX` repository by typing: ::

    git clone https://github.com/AMReX-Codes/amrex.git

This will create a folder called ``amrex/`` on your machine. Set the environment variable, ``AMREX_HOME``, on your
machine to point to the path name where you have put `AMReX`::

        export AMREX_HOME=/path/to/amrex/
        
3. Clone the `Pele` repository: ::

    git clone git@github.com:AMReX-Combustion/PelePhysics.git

This will create a folder called ``PelePhysics`` on your machine.
Set the environment variables, ``PELE_PHYSICS_HOME``, to where you put these.

4. Periodically update both of these repositories by typing ``git pull`` within each repository.

