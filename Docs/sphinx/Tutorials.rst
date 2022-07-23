.. role:: cpp(code)
   :language: c++

.. _sec:tutorials:

Tutorials
=========

This sections includes several tutorials - TODO

.. toctree::
   :maxdepth: 1

Tutorial 1 - Generating NC12H26 QSS mechanism with analytical jacobian
----------------------------------------------------------------------

Update ``poetry``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry update  

Make sure the list of non-QSS species is correct in ``Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml``

Generate ``qssa.yaml``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml

Generate ``mechanism.H`` and ``mechanism.cpp``::
  
  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa.yaml --format_input qssa_input.toml --symbolic_jacobian

Tutorial 2 - Generating NC12H26 QSS mechanism without analytical jacobian
-------------------------------------------------------------------------

Update ``poetry``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry update  

Make sure the list of non-QSS species is correct in ``Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml``

Generate ``qssa.yaml``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml

Generate ``mechanism.H`` and ``mechanism.cpp``::
  
  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa.yaml

