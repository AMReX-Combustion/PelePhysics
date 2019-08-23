.. highlight:: rst

`PelePhysics` Quickstart
============================

Greetings impatient user. Once again, note that this documentation focuses on the CHEMISTRY part of `PelePhysics`.

- If you are already familiar with this code, have it installed already and would simply like to know which chemistry-related keywords and/or environment variables to set in your various input files to perform a simulation with one of the codes available in the `PeleSuite`, then I invite you to directly skip to section :ref:`sec:subsPPOptions`. 
- If you are a complete beginner, I urge you to carefully read the two following chapters :ref:`sec:GetPP` and :ref:`sec:GetCVODE`, to properly set-up your working environment.
- If you are in a hurry but still would like more context, visit section :ref:`sec:subsWD` to be referred to portions of this document that are of interest to you.


.. _sec:GetPP:

Obtaining `PelePhysics`
----------------------------


First, make sure that "Git" is installed on your machine---we recommend version 1.7.x or higher. Then...

1. Download the `AMReX` repository by typing: ::

    git clone https://github.com/AMReX-Codes/amrex.git

This will create a folder called ``amrex/`` on your machine. Set the environment variable, ``AMREX_HOME``, on your machine to point to the path name where you have put `AMReX`::

        export AMREX_HOME=/path/to/amrex/
        
2. Clone the `Pele` repository: ::

    git clone git@github.com:AMReX-Combustion/PelePhysics.git

This will create a folder called ``PelePhysics`` on your machine. Set the environment variable ``PELE_PHYSICS_HOME`` to where you put these.

3. Periodically update both of these repositories by typing ``git pull`` within each repository.


.. _sec:GetCVODE:

Install CVODE and `SuiteSparse`
---------------------------------------
**The user is in charge of installing the proper CVODE version, as well as installing and properly linking the KLU library if sparsity features are needed.**


SuiteSparse
^^^^^^^^^^^^^^^^^^^^^^^^
`SuiteSparse` is a suite of Sparse Matrix software that is very handy when dealing with big kinetic mechanisms (which are usually very sparse). 
In such a case, CVODE can make use of the KLU library, which is part of `SuiteSparse`, to perform sparse linear algebra.
Documentation and further information can be found `here <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_. 

At the time this note is written, the recommended **SuiteSparse version** is **5.4.0**. Follow these steps to set-up your working environment and build the required libraries:

1. Go to `the website <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_ and download the compressed file for the recommended version
2. Copy the tar file into ``$PELE_PHYSICS_HOME/ThirdParty``
3. Untar ('tar -zxvf'), cd into it and type 'make' into the following folders: ``SuiteSparse_config``, ``AMD``, ``COLAMD``, ``BTF``
4. Go into ``metis-5.1.0`` and type 'make config shared=1' followed by 'make'
5. Go into ``KLU`` and make
6. Check that all dynamic libraries have correctly been generated and copied into the folder ``$PELE_PHYSICS_HOME/ThirdParty/SuiteSparse/lib`` 
7. It is recommended that you add the path ``$PELE_PHYSICS_HOME/ThirdParty/SuiteSparse/lib`` to your ``LD_LIBRARY_PATH``, for precaution
8. Note that depending upon your compiler, the static ``.a`` versions of the libraries might also be required. In such a case, you can copy them from each folder into ``lib``

CVODE
^^^^^^^^^^^^^^^^^^^^^^^^
CVODE is a solver for stiff and nonstiff ordinary differential equation (ODE) systems. Documentation and further information can be found `online <https://computing.llnl.gov/projects/sundials/cvode>`_.
At the time this note is written, the recommended **CVODE version** is **v5.0.0-dev.1**. 

The CVODE sources are distributed as compressed archives, with names following the convention ``cvode-x.y.z.tar.gz``. They can be downloaded by following 
`this link <https://computation.llnl.gov/projects/sundials/sundials-software>`_.  However, we have designed a script enabling to install the correct version the correct way. Simply:

1. Go into ``$PELE_PHYSICS_HOME/ThirdParty`` 
2. Execute either ``get_sundials_v5dev1.sh`` or ``get_sundials_v5dev1_CUDA.sh`` depending on your application and machine
3. Set the ``CVODE_LIB_DIR`` environment variable to point towards the location where all CVODE libraries have been generated. If you followed these guidelines, it should be ``$PELE_PHYSICS_HOME/ThirdParty/sundials/instdir/lib/`` 

Note that if you do not want to use the KLU library, you can also disable the flags (``-DKLU_ENABLE``) in the scripts. 

