.. highlight:: rst

.. role:: cpp(code)
   :language: c++

CVODE implementation in `PelePhysics`
======================================

.. _sec:subsDiffReacts:

The different reactors
-----------------------

Throughout this document, what we call a `reactor` is in fact a zero-dimensional model, 
representing the simplest form of a chemically reacting system. Depending upon the choice of state variables 
driving the system, several different types of reactor can be considered; 
and the "correct" choice is case dependent. In general, the state variables for a reactor model are

- The reactor mass
- The reactor volume
- The energy of the system
- The mass fractions for each species

The most common type of reactor is the `constant-volume` (CV) reactor, which is the one used to advance the chemistry 
within `PeleC`. This reactor type is equivalent to a rigid vessel with fixed volume but variable pressure. 
In `PelePhysics`, the constant-volume constraint is ensured by keeping the density :math:`\rho` fixed 
-since there is no change of mass; and the indirect choice of energy in the CV reactor implementation is the total energy 
:math:`E`. :math:`E`'s evolution in our case is solely due to a constant external source term :math:`\dot{E}_{ext}`, which accounts 
for the effects of advection and convection in the Spectral Deferred Correction (SDC) scheme that all Pele codes use. 
In that sense, the CV reactor is an abstraction and is not a true closed vessel.

Note that CVODE still integrates the mass fractions (:math:`\rho Y`) together with energy for stability reasons, 
but a change of variable is applied to effectively transport the temperature :math:`T` via

.. math::

    \rho C_v \frac{\partial T}{\partial t} = \rho\dot{E}_{ext}  - \sum_k e_k {\dot{\omega}_k}^M

where the :math:`e_k` are the species internal energy and :math:`{\dot{\omega}_k}^M` is the species :math:`k` mass production rate. 

In a second implementation, that we will label `constant-volume-enthalpy` (CVH), the mass-weighted total enthalpy :math:`\rho H` is used and 
conserved along with :math:`\rho`. This reactor type is also an abstraction. Here also, :math:`\rho H` 
evolves according to an external source term :math:`\dot{\rho H}_{ext}`, and in CVODE, the mass fractions (:math:`\rho Y`) and 
temperature :math:`T` are integrated according to

.. math::

    \rho C_p \frac{\partial T}{\partial t} = \rho\dot{H}_{ext}  - \sum_k h_k  {\dot{\omega}_k}^M

where the :math:`h_k` are the species internal energy. 

.. _sec:subsubValidCVreact:

Validation of the CV reactor implementation (with CANTERA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`CANTERA <https://cantera.org/>`_ is an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes. 
It is a very robust and fast tool written in C++ that is also based on CVODE to perform the chemistry integration. 
CANTERA is well recognized in the combustion community, and by comparing our results to CANTERA reactor simulations, 
we will be able to validate our implementation. 

Note that only the CV reactor model described above can be validated, since as we discussed before, 
the CVH reactor model is an abstraction needed for our Low-Mach PeleLM chemistry integration. Also, to have a real CV reactor, 
the external source terms on the energy and species in `PelePhysics` have been set to 0 (see :ref:`sec:subsDiffReacts`).

The parameters chosen to initialize the simulation in both CANTEA and `PelePhysics` are described in 
Table :numref:`tab:ReactEvalCVODE`. Two sets of runs are performed, with two mechanisms available in `PelePhysics`. 
Note that small sub-steps are taken until the final time is reached, 
but CVODE's internal machinery can subdivides the :math:`dt` even further. 
For the purpose of validation, the direct dense solver of CVODE is selected 
in `PelePhysics` (see section :ref:`sec:subsPPOptions`).

.. _tab:ReactEvalCVODE:

.. table::
    :align: center

    +------------+-----------------+-------------+----------------+-------------+----------------+-----------------+
    | Mechanism  |     Mixture     |  Initial T  |  Initial phi   |   Pressure  |       dt       |    Final time   |
    +------------+-----------------+-------------+----------------+-------------+----------------+-----------------+
    |  Li Dryer  |      H2/O2      |   1500 K    |      0.8       |  101325 Pa  |     1.0e-8s    |     1.0e-5s     |
    +------------+-----------------+-------------+----------------+-------------+----------------+-----------------+

Results are plotted in Fig :numref:`fig:ReactEvalCVODE` and :numref:`fig:ReactEvalCVODESpecs`. for the :math:`H_2/O_2` mixture. 
All curves are indistinguishable, so the relative error of all major quantities is also plotted in Fig. :numref:`fig:ReactEvalCVODEErrss`. 
Note that :math:`H_2` and :math:`O_2` relative errors have similar features, and that relative errors observed 
for :math:`H` and :math:`H_2O` are representative of those exhibited by, respectively, intermediates and products.

.. |a| image:: ./Visualization/Main.001.png
     :width: 100%

.. |b| image:: ./Visualization/Specs.png
     :width: 100%

.. |c| image:: ./Visualization/ERRs.png
     :width: 100%

.. _fig:ReactEvalCVODE:

.. table:: Evolution of temperature, pressure and enthalpy in a CV reactor, computed with the LiDryer mechanism. Black: CANTERA, red: PelePhysics.
     :align: center

     +-----+
     | |a| |
     +-----+
..
    .. figure:: ./Visualization/Main.001.png
     :width: 100%
     :name: fig-ReactEvalCVODE
     :alt: Evolution of temperature, pressure and enthalpy in a CV reactor, computed with the LiDryer mechanism. Black: CANTERA, red: PelePhysics.


.. _fig:ReactEvalCVODESpecs:

.. table:: Evolution of major species in a CV reactor, computed with the LiDryer mechanism. Black: CANTERA, red: PelePhysics. 
     :align: center

     +-----+
     | |b| |
     +-----+

.. _fig:ReactEvalCVODEErrss:

.. table:: Relative errors on the temperature, pressure, enthalpy and major species in a CV reactor, computed with the LiDryer mechanism. 
     :align: center

     +-----+
     | |c| |
     +-----+


Overall, considering the many CVODE controlling parameters, results are deemed acceptable and that 
concludes the validation of the reactors implementations in `PelePhysics`.


.. _sec:subsPPOptions:

Activating the different options via the input files
-----------------------------------------------------

Choosing between DVODE/CVODE is done at compile time, 
via the ``GNUmakefile``. On the other hand, the type of reactor and numerical algorithm 
are selected via keywords in the input file. There is a subtlety though: 
when any sparsity feature is required, the choice should also be made a compile time; 
and if the option is not selected, subsequent options via keywords can either lead to an error or fall back to a dense formulation 
of the problem. This is discussed in more depth in what follows.

Preliminary step: CVODE and SuiteSparse
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**The user is in charge of installing the proper CVODE version, as well as installing and properly linking the KLU library if sparsity features are needed.**

As previously said, the **CVODE version** employed to perform all tests presented in this document is **v4.0.2**. 
The CVODE sources are distributed as compressed archives, with names following the convention ``cvode-x.y.z.tar.gz``. They can be downloaded by following 
`this link <https://computation.llnl.gov/projects/sundials/sundials-software>`_. 
The first step is to extract the files: ::

    % tar xzf cvode-x.y.z.tar.gz

This will extract source files under a directory ``cvode-x.y.z``. An ``INSTALL_GUIDE.pdf`` will be present in these source files, 
explaining how to properly install and build the needed librairies. We strongly recommend following the convention and generating 2 directories: :: 

    % cd cvode-x.y.z
    % mkdir builddir instdir
    % cd builddir

Now, the install directory defaults to ``/usr/local/instdir`` but should be changed to point to the ``cvode-x.y.z/instdir`` folder you just generated by setting the
``-DCMAKE_INSTALL_PREFIX`` flag (this is all explained in the ``INSTALL_GUIDE.pdf``)

If the KLU library is to be used, then the proper flags should be set. With the version sub-mentioned, those are 
``-DKLU_ENABLE`` (should be set to ``ON``)  ``-DKLU_INCLUDE_DIR`` (should point to the location where your ``klu.h`` file is) and ``-DKLU_LIBRARY_DIR`` (which should point to the location where your KLU library was generated). 

The KLU library is part of the `SuiteSparse` distribution that can be downloaded `here <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_.
The SuiteSparse version used for all tests presented in this document is **v4.5.0**, and we do not guarantee that the linking will properly work with subsequent releases. 
A ``README`` file will explain how to generate the KLU library.

Once CVODE has been build and installed -with proper linking to the KLU library if necessary, 
the ``CVODE_LIB_DIR`` environment variable should be set, pointing towards the location where all CVODE libraries have been generated.
If the convention has been respected, this should be ::

    CVODE_LIB_DIR=pathToYourcvode-x.y.z/instdir/lib/

.. _subsubs:GNUtype:

The GNUmakefile
^^^^^^^^^^^^^^^^^^^^^^^^

The default setting is to use DVODE in `PelePhysics`; i.e, if no modifications are done to the original ``GNUmakefile``, 
then this option should automatically be selected. To activate CVODE, the following line should be added: ::

    USE_SUNDIALS_3x4x = TRUE

Note that this is an AMREX flag, so it will automatically be recognized throughout the `Pele` codes. 

If sparsity features are required, then the following line should also be added:

.. code-block:: c++

    USE_KLU = TRUE
    ifeq ($(USE_KLU), TRUE)  
        DEFINES  += -DUSE_KLU
    endif

In this case, the location of the ``klu.h`` should be specified via: ::

    INCLUDE_LOCATIONS += path_to_your_klu.h

The only place where this flag will be used is in the C++ file where all CVODE calls are done.


The input file
^^^^^^^^^^^^^^^^^^^^^^^^

If DVODE has been enabled in the ``GNUmakefile``, no modifications of the input file are required. If CVODE is selected, 
then the original input file will trigger a dense direct solve without Analytical Jacobian. 
Three keywords control the algorithm.

- ``cvode_iE`` enable to switch from a CV reactor (``=1``) to a CVH reactor (``=2``).
- ``ns.cvode_iDense`` controls the numerical method: choose ``1`` to enable the dense direct linear solver, 
  ``5`` for the sparse direct linear solver (if the KLU library has been linked) and ``99`` for the Krylov iterative solver
- ``ns.cvode_iJac`` is a bit less obvious. If ``ns.cvode_iDense = 1``, then ``ns.cvode_iJac = 1`` will activate 
  the use of an Analytical Jacobian. But if ``ns.cvode_iDense = 99``, then ``ns.cvode_iJac = 1`` will activate 
  the preconditioned GMRES solver while ``ns.cvode_iJac = 0`` will activate the non-preconditioned GMRES solver. 
  Additionally, if ``ns.cvode_iDense = 99``, ``ns.cvode_iJac = 1`` and the KLU library is linked, then the preconditioned solve 
  is done in a sparse format. Note that with ``ns.cvode_iDense = 5``, the only allowed option is ``ns.cvode_iJac = 1``.


.. _sec:subsReactEvalCvode:

The ReactEvalCVODE test case in details
---------------------------------------

This tutorial has been adapted from the ReactEval tutorial employed in the series of regression tests to monitor the DVODE chemistry integration. 
The domain considered is a :math:`2x1024x2` box, where the initial temperature is different in each :math:`(i,j,k)-` cell, according to a :math:`y-` evolving sinusoidal profile, see Fig. :numref:`fig:ErrH2`:

.. math::

    T(i,j,k) =  T_l + (T_h-T_l)\frac{y(i,j,k)}{L} + dTsin\left(2\pi\frac{y(i,j,k)}{P}\right) 

The different parameters involved are summarized in Table :numref:`tab::ParamReactEvalCvode`. The initial composition 
is the same in every cell, and is a mixture of 0.1 :math:`CH_4`, 0.2 :math:`O_2` and 0.7 :math:`N_2` in mass fractions. 
The initial pressure is 1 atm and the mechanism considered is the DRM mechanism.

.. _tab::ParamReactEvalCvode:

.. table:: Parameters used to initialize T in the ReactEvalCVODE test case
    :align: center

    +------------+-----------------+-------------+----------------+-------------+
    | Tl         |     Th          |  dT         |  L             |   P         |
    +------------+-----------------+-------------+----------------+-------------+
    |  2000 K    |      2500 K     |   100 K     |      1024      |  L/4        |
    +------------+-----------------+-------------+----------------+-------------+


.. _fig:ErrH2:

.. figure:: ./Visualization/Case_ReactEvalCvode.001.png
     :width: 50%
     :align: center
     :name: fig-ReactEvalCVODE
     :target: ./Visualization/Case_ReactEvalCvode.001.png
     :alt: The ReactEvalCVODE test case

     The ReactEvalCVODE test case

The GNUmakefile
^^^^^^^^^^^^^^^^^^^^^^^^

For this example, the ``USE_SUNDIALS_3x4x`` flag should be set to true, as the ODE integration 
is called from the C++ routine directly and will thus not work with DVODE (only accessible via FORTRAN with the current implementation).
Additionally, the ``FUEGO_GAS`` flag should be set to true and the chemistry model should be set to ``drm19``. The full file reads as follows:

.. code-block:: c++

    PRECISION  = DOUBLE                                                                                                                   
    PROFILE    = FALSE
    
    DEBUG      = TRUE
    
    DIM        = 3
    
    COMP       = gcc
    FCOMP      = gfortran
    
    USE_MPI    = TRUE
    USE_OMP    = FALSE
    
    FUEGO_GAS  = TRUE
    
    # define the location of the PELE_PHYSICS top directory
    PELE_PHYSICS_HOME    := ../../../..
    
    #######################
    USE_SUNDIALS_3x4x = TRUE
    
    USE_KLU = FALSE
    ifeq ($(USE_KLU), TRUE)
        DEFINES  += -DUSE_KLU
        include Make.CVODE
    endif
    
    #######################
    ifeq ($(FUEGO_GAS), TRUE)
      Eos_dir         = Fuego
      Chemistry_Model = drm19
      Reactions_dir   = Fuego
      Transport_dir   = Simple
    else
      Eos_dir       = GammaLaw
      Reactions_dir = Null
      Transport_dir = Constant
    endif
    
    Bpack   := ./Make.package
    Blocs   := .

    include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics         

where the ``Make.CVODE`` contains the link to the `SuiteSparse` include files (see :ref:`subsubs:GNUtype`).


The input file
^^^^^^^^^^^^^^^^^^^^^^^^

The run parameters that can be controlled via the input files are as follows: ::

    dt = 1.e-05  
    ndt = 100
    cvode_ncells = 1
    
    cvode_iE = 1
    ns.cvode_iDense = 1
    ns.cvode_iJac = 0

    amr.plot_file       = plt

so in this example, a **CV reactor model is chosen** to integrate each cell, and the **dense direct solve without analytical Jacobian** is activated. 
Each cell is then integrated for a total of :math:`1.e-05` seconds, with 100 external time steps. 
This means that the actual :math:`dt` is :math:`1.e-07s`, a little more than what is used in the `PeleC` code, 
but consistent with what will be used in `PeleLM`. The number of cells to be integrated simultaneously is 1 [#Foot1]_.


Results
^^^^^^^^^^^^^^^^^^^^^^^^

It took 2m44s to integrate the 4096 cells of this box, with 4 MPI processes and no OMP process. 
The resulting temperature evolution for all cells is displayed in Fig. :numref:`fig:ReacEvalCv`.


.. _fig:ReacEvalCv:

.. figure:: ./Visualization/ReactEvalCv.001.png
     :width: 100%
     :align: center
     :name: fig-ReactEvalCv
     :alt: Evolution of temperature in the 2x1024x2 example box, using a CV reactor and a dense direct solve, and computed with the DRM mechanism. Black: $t=0$, red: $t=1e-05s$

     Evolution of temperature in the 2x1024x2 example box, using a CV reactor and a dense direct solve, and computed with the DRM mechanism. Black: t=0s, red: t=1e-05s


To go further: ReactEvalCVODE with the KLU library
---------------------------------------------------

The GNUmakefile
^^^^^^^^^^^^^^^^^^^^^^^^

Only the middle part of the ``GNUmakefile`` is modified compared to the previous example, section :ref:`sec:subsReactEvalCvode`.

.. code-block:: c++

    ...
    #######################
    USE_SUNDIALS_3x4x = TRUE
    
    USE_KLU = TRUE
    ifeq ($(USE_KLU), TRUE)
        DEFINES  += -DUSE_KLU
        include Make.CVODE
    endif
    
    #######################
    ...


The input file
^^^^^^^^^^^^^^^^^^^^^^^^

For the KLU library to be of use, a solver utilizing sparsity features should 
be selected. We modify the input file as follows: ::

    dt = 1.e-05  
    ndt = 100
    cvode_ncells = 1
    
    cvode_iE = 1
    ns.cvode_iDense = 99
    ns.cvode_iJac = 1
    
    amr.plot_file       = plt

So that now, a preconditioned iterative Krylov solver is selected, where the preconditioner is specified in a sparse format.

Results
^^^^^^^^^^^^^^^^^^^^^^^^

This run now takes 3m29s to run. As expected from the dense Jacobian of the system obtained when using the small DRM mechanism 
(the fill in pattern is :math:`>90 %`), using an iterative solver does not enable to reach speed-ups over the simple dense direct 
solve. **NOTE**, and this is important, that this tendency will revert when sufficiently small time steps are used. 
For example, if instead of :math:`1e-7s` we took time steps of :math:`1e-8s` (consistent with `PeleC` time steps), then using 
the iterative GMRES solver would have provided significant time savings. This is because the smaller the time step the 
closer the system matrix is from the identity matrix and the GMRES iterations become really easy to complete.

This example illustrates that choosing the "best" and "most efficient" algorithm is far from being a trivial task, 
and depends upon many factors. Table :numref:`tab:RunsReactEvalCvode` provides a summary of the run time in solving the 
ReactEvalCVODE example with the various available CVODE linear solvers.

.. _tab:RunsReactEvalCvode:

.. table:: Summary of ReactEvalCvode runs with various algorithms
    :align: center

    +------------+-----------------+-------------+----------------+----------------+-----------------+
    |  Solver    |     Direct      |  Direct     |  Direct        |   Iter.        |   Iter.         |
    |            |     Dense       |  Sparse AJ  |  Dense AJ      |   Precond. (D) |   Precond. (S)  |
    +------------+-----------------+-------------+----------------+----------------+-----------------+
    |  KLU       |     OFF         |  ON         |  OFF           |   OFF          |   ON            |
    +============+=================+=============+================+================+=================+
    |            |     Dense       |   Sparse AJ |     Dense AJ   |   Precond. (D) |   Precond. (S)  |
    +------------+-----------------+-------------+----------------+----------------+-----------------+

..
    \hline 
    "cvode_iE"              & 1                      &     1                     & 1                                    & 1                        &  1      \\
    "ns.cvode_iDense" & 1                      &    5                      & 1                                    & 99                      &  99 \\ 
    "ns.cvode_iJac"      & 0                      &    1                      & 1                                    & 1                        &  1 \\ 
    \hline
    Run time $[s]$      &  2m44s           & 1m59s                      & 1m42s                           & 3m05s                    &  3m29s   \\
    \end{tabular}
    


Current Limitations
--------------------



Tricks and hacks, stuff to know
---------------------------------



How does CVODE compare with DVODE ?
-----------------------------------

.. [#Foot1] NOTE that only one cell at a time should be integrated with CVODE right now. The vectorized version on CPU is still WIP and not properly implemented for all linear solver.
