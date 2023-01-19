.. role:: cpp(code)
   :language: c++

.. _sec:tutorials:

Tutorials
=========

This sections includes several tutorials - TODO

.. toctree::
   :maxdepth: 1

Tutorial 1 - Generating NC12H26 QSS mechanism with analytical Jacobian
----------------------------------------------------------------------
.. _sec_tutqss1:

Update ``poetry``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry update  

Make sure the list of non-QSS species is correct in ``Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml``

Generate ``qssa.yaml``::

  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml

Generate ``mechanism.H`` and ``mechanism.cpp``::
  
  $ cd ${PELE_PHYSICS_HOME}/Support/ceptr
  $ poetry run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa.yaml --qss_format_input ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa_input.toml --qss_symbolic_jacobian

We recommend using the following ``qssa_input.toml``::
 
  $ [Readability]
    hformat =                       "gpu"
    
    [Arithmetic]
    remove_1 =                      true
    remove_pow =                    true
    remove_pow10 =                  true
    
    [Replacement]
    min_op_count =                  0
    min_op_count_all =              10
    gradual_op_count =              true
    remove_single_symbols_cse =     true
    
    [Recycle]
    store_in_jacobian =             true
    recycle_cse =                   true
    
    [Characters]
    round_decimals =                true


Tutorial 2 - Generating NC12H26 QSS mechanism without analytical Jacobian
-------------------------------------------------------------------------
.. _sec_tutqss2:

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

