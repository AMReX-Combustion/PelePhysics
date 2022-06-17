# TO DO
I tried to not modify functions that did not need be modified but we should probably revamp this to make the code more readable. For ex: functions like `sorted_kc_exp` now throw back 1 or 2 output if symbolic recording is enabled or not. We should clean that

1. Homogeneize functions (see above)
2. ~~Propagate differentiation through `comp_sc_qss`~~
3. ~~Verify that `qf_qss` and `qr_qss` are consistent with `comp_qss_coeff`~~
4. ~~Verify that `sc_qss` is consistent with `comp_sc_qss`~~ See `debugDodecaneLuQss`
5. ~~Make sure that the expressions printed are compact~~
6. ~~Make sure code is correctly formatted~~
7. Verify that `dsc_qss/dsc` is correct compared to finite difference
8. Modify `aJacobian` to include the dependence of `sc_qss` with `sc` and `T`.
  - ~~Symbolic encoding of `wdot`~~
  - ~~Verify that `wdot` values obtained with `sympy`~~
  - Chain rule to get `dscqss_dsc` -> Ongoing
  - Verify `dscqss_dsc` -> Ongoing (some `dscqss_dsc` are correct but others arent)
  - Chain rule to get `dwdot_dsc`
  - ~~Finite difference to get `dwdot_dT`~~
  - Verify `dTdot` terms -> Ongoing
9. Ensure computational efficiency
  - ~~Use `symEngine` when computing `diff`~~
  - Make sure we do not recompute twice the same `sme.diff` -> Ongoing
10. Try 0D example (`ReactEval_dodecanelu_qss`) with and without Jacobian fix

# PelePhysics
*A respository of physics databases and implementation code for use with the `Pele` suite of of codes*

## CEPTR

This `PelePhysics` repository contains the CEPTR source code generation tool in order to support the inclusion of chemical models specified in the Cantera yaml format. CETPR derives from FUEGO, which was originally created by Michael Aivazis at CalTech, and donated to CCSE in 2001.  Originally, FUEGO was part of a larger Python-based set of workflow tools Michael developed for managing and publishing simulation-based studies.  FUEGO itself was developed as a drop-in replacement for the CHEMKIN library, and provided "hand-coded" replacement routines for evaluation of EOS thermodynamic functions that are considerably more efficient than their CHEMKIN counterparts.  Since 2001, CCSE has continued to modify FUEGO independently for its own use so that the current version here bears little resemblance to Michael's original code, or its intentions. CEPTR has adapted FUEGO code to use the Cantera yaml format and is now the preferred way of generating mechanism files.

Typically, Cantera *mechanisms* (combustion models) are transmitted via a set of yaml files, which provide data for
* *Species* - definition of chemical species, as composed by fundamental elements
* *Chemistry* - modified Arrhenius reaction parameters, including third body efficiencies and pressure dependencies
* *Thermodynamics* - NASA-formatted polynomial expressions, assuming a mixture of ideal gases
* *Transport* - physical parameters for each of the chemical species required for the `EGLib` (or comparable) library for computing species- and temperature-dependent molecular transport coefficients

This information is read by the CEPTR tool, and source code files are generated that are to linked directly into an executable in the `Pele` code suite in order to provide data and functions necessary to evaluate chemistry, transport and thermodynamic relationships. `PelePhysics` supports the following in addition to standard Cantera style models:
* *Soave-Redlich-Kwong* real gas equation of state, implemented to be consistent with corresponding ideal gas model at suitable conditions
* *Gamma-law* gas equation of state
* *Constant coefficient* transport
* *Simple* A simplified mixture-averaged transport model consistent with EGLib under suitable conditions
* *Inert* Chemistry setting to ignore reactions
* Utilization of auxiliary source files (C++) that follow the Cantera species production rate function interface (e.g., QSS approximations)

If the `Pele` codes are built with `Eos_dir = Fuego`, the make system variable `Chemistry_Model` must be set to one of the models (subfolders) that exist in the `${PELE_PHYSICS_HOME}/Support/Mechanism/Models` folder. The repository currently provides (multiple) models for each of the following mixtures:
* air
* hydrogen
* methane
* ethane
* *n*-dodecane
* di-methyl ether
* ... more can be added (see below)


### NOTE: Non-ideal equations of state

`PelePhysics` currently supports a cubic EOS model: `Soave-Redlich-Kwong`.  It is built on top of the ideal gas models, and is selected by specifying its name as the `Eos_dir` during the build (the make system requires that both `Eos_dir` and `Chemistry_Model` be specified).  Any additional parameters required for the EOS (e.g., attractions, repulsions, critical states) are either included in the underlying CEPTR database used to generate the source file model implementation, or else are inferred from the input model data.

## Model generation procedures
This repository provides the tools necessary for generating new Pele-compatible combustion mechanisms. Please refer to the `CEPTR documentation <https://amrex-combustion.github.io/PelePhysics/Ceptr.html>`_ for instructions on generating mechanism models. Make sure that you edit the `GNUmakefile` where you want to use this (in, e.g., `PeleC/Exec`) so that `Chemistry_Model` is `XXX`.  In `PeleC/Exec/Make.PeleC`, the model is expected to be in the folder `${PELE_PHYSICS_HOME}/Support/Mechanism/Models/$(Chemistry_Model)`, and it is expected that the folder contains a `Make.package` file to include, so make sure things are where they need to be. Refer to other mechanisms for additional guidance.
