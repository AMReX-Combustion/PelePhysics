.. highlight:: rst

.. _SootInputs:

Soot Flags and Inputs
======================

* In the ``GNUmakefile``, specify ``USE_SOOT = TRUE`` and ``NUM_SOOT_MOMENTS = N`` where ``N`` is the number of moments to be used in the solution, either 3 or 6.

* Depending on the gas phase solver, soot solving functionality can be turned on in the input file using ``pelec.add_soot_src = 1`` or ``peleLM.do_soot_solve = 1``.

* The chemistry model specified with ``Chemistry_model =`` in  ``GNUmakefile`` must contain the PAH inception species, as well as ``H2``, ``H``, ``OH``, ``H2O``, ``CO``, ``C2H2``, and ``O2``.

* A PAH inception species must be provided in the input file using ``soot.incept_pah =``. Currently, only one of three inputs are accepted: ``A2``, ``A3``, or ``A4``; which correspond to naphthalene (C10H8), phenathrene (C14H10), or pyrene (C16H10).

  * If the inception species is named something other than ``A#`` in the chemistry model, a different name can be specified using ``soot.pah_name =``. However, ``soot.incept_pah`` must be set to ``A2``, ``A3``, or ``A4``.


