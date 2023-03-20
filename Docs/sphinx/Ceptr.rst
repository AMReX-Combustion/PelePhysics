.. highlight:: rst

.. _sec:ceptr:
               
CEPTR: Chemistry Evaluation for Pele Through Recasting
======================================================

We use CEPTR to generate C++ mechanism files from `Cantera <https://cantera.org>`_ yaml chemistry files. CEPTR is a python package part of the PelePhysics source code.

Software requirements
---------------------

The CEPTR package uses `poetry <https://python-poetry.org/docs/#installation>`_ to manage the Python dependencies. Poetry is therefore required to use CEPTR and can typically be installed through system package managers (e.g. HomeBrew) or following the instructions in poetry's documentation.

To install CEPTR dependencies::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry update


Usage
-----

There are three ways to use CEPTR to generate C++ mechanism files for a given chemistry

1. Using CEPTR directly::

     $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
     $ poetry run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/LiDryer/mechanism.yaml

2. Using a helper script in the directory containing the ``mechanism.yaml`` file::

     $ ./convert.sh

3. Using a helper script in the ``Models`` directory::

     $ bash ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/converter.sh -f ./LiDryer/mechanism.yaml


For non-reduced chemistries, CEPTR can take a file with a list of ``mechanism.yaml`` files to convert::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run convert -l ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/list_mech

For reduced chemistries, CEPTR can take a file with a list of ``qssa.yaml`` and ``qssa_input.toml`` to convert::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run convert -lq ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/list_qss_mech

For generating ``qssa.yaml`` for reduced chemistries, CEPTR can take a file with a list of ``skeletal.yaml`` and ``non_qssa_list.yaml``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run qssa -lq ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/list_qss_mech


Converting CHEMKIN files
------------------------
.. _sec_convertCK:

We rely on Cantera's ``ck2yaml`` utility to convert CHEMKIN files to the Cantera yaml format (and proceed as above with CEPTR on the resulting yaml file)::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run ck2yaml --input=${PATH_TO_CHEMKIN_DIR}/mechanism.inp --thermo=${PATH_TO_CHEMKIN_DIR}/therm.dat --transport=${PATH_TO_CHEMKIN_DIR}/tran.dat --permissive

The files ``tran.dat`` and ``therm.dat`` are optional if already included in the ``.inp`` file.

Generating a QSS chemistry file
-------------------------------

To generate a QSS chemistry yaml file from another yaml file, one executes::

  $ poetry run qssa -f ${PATH_TO_YAML}/skeletal.yaml -n ${PATH_TO_YAML}/non_qssa_list.yaml

The full list of options is::

  $ poetry run qssa -h
  usage: qssa [-h] -f FNAME -n NQSSA [-m {0,1,2}] [-v]

  Mechanism converter

  optional arguments:
    -h, --help            show this help message and exit
    -f FNAME, --fname FNAME
                          Mechanism file
    -n NQSSA, --nqssa NQSSA
                          Non-QSSA species list
    -m {0,1,2}, --method {0,1,2}
                          QSSA method (default: 2)
    -v, --visualize       Visualize quadratic coupling and QSSA dependencies

For a detailed description of these options and a further information on the way QSS mechanism are treated in `CEPTR` the reader may consult :ref:`the QSS section <sec_qss>`.

See Tutorials (:ref:`Generating NC12H26 QSS mechanism with analytical jacobian <sec_tutqss1>` and :ref:`Generating NC12H26 QSS mechanism without analytical jacobian <sec_tutqss2>`) for generating QSS mechanisms from the ``.yaml`` files.


