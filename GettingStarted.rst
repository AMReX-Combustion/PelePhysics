.. highlight:: rst

.. _GettingStarted:

Getting Started
===============

The `PeleMP` code can be obtained via `GitHub <https://github.com/AMReX-Combustion/PeleMP>`_.
`PeleMP` is the repository of multiphysics capabilities for use with other `Pele` codes:

* `PeleC <https://amrex-combustion.github.io/PeleC/>`_
* `PeleLM <https://amrex-combustion.github.io/PeleLM/>`_
* `PeleLMeX <https://amrex-combustion.github.io/PeleLMeX/>`_


Since the suite of `Pele` codes are in active development, the easiest method to access `PeleMP` for end-users is through the `PeleProduction` repository.

1. Clone the `PeleProduction` repo and switch the appropriate branch: ::

     git clone https://github.com/AMReX-Combustion/PeleProduction.git PeleProduction
     cd PeleProduction
     git checkout Landon/soot_and_spray
     git submodule init
     git submodule update

2. Access the directory of test problems for `PeleMP` ::

     cd PeleMPruns

3. Access the directory of the `Pele` code you wish to interface with, either `PeleC`, `PeleLM`, or `PeleLMeX`.

4. Access the directory of the relevant problem type, either ``spray_flame`` for a spray injection case, or ``laminar_soot`` for a sooting flame case. These are just two example problems to help demonstrate the capabilities of `PeleMP`.

5. Modify the ``GNUmakefile`` according to the desired build type (compiler, MPI, chemistry, soot and spray flags, etc).

6. Run the build script in the directory ``./build.sh``.

.. toctree::
   :maxdepth: 1

   SprayInputs


