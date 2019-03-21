# PelePhysics
*A respository of physics databases and implementation code for use with the `Pele` suite of of codes*

## FUEGO

This `PelePhysics` repository contains the FUEGO source code generation tool in order to support the inclusion of chemical models specified in the CHEMKIN-II format.  Fuego was originally created by Michael Aivazis at CalTech, and donated to CCSE in 2001.  Originally, FUEGO was part of a larger Python-based set of workflow tools Michael developed for managing and publishing simulation-based studies.  FUEGO itself was developed as a drop-in replacement for the CHEMKIN library, and provided "hand-coded" replacement routines for evaluation of EOS thermodynamic functions that are considerably more efficient than their CHEMKIN counterparts.  Since 2001, CCSE has continued to modify FUEGO independently for its own use so that the current version here bears little resemblance to Michael's original code, or its intentions.

Typically, CHEMKIN-based *mechanisms* (combustion models) are transmitted via a set of ASCII files, which provide data for
* *Species* - definition of chemical species, as composed by fundamental elements
* *Chemistry* - modified Arrhenius reaction parameters, including third body efficiencies and pressure dependencies
* *Thermodynamics* - NASA-formatted polynomial expressions, assuming a mixture of ideal gases
* *Transport* - physical parameters for each of the chemical species required for the `EGLib` (or comparable) library for computing species- and temperature-dependent molecular transport coefficients

This information is read by the FUEGO tool, and source code files are generated that are to linked directly into an executable in the `Pele` code suite in order to provide data and functions necessary to evaluate chemistry, transport and thermodynamic relationships. `PelePhysics` supports the following in addition to standard CHEMKIN style models:
* *Soave-Redlich-Kwong* real gas equation of state, implemented to be consistent with corresponding ideal gas model at suitable conditions
* *Gamma-law* gas equation of state
* *Constant coefficient* transport
* *Simple* A simplified mixture-averaged transport model consistent with EGLib under suitable conditions
* *Inert* Chemistry setting to ignore reactions
* Utilization of auxiliary source files (C++,f90,c) that follow the CHEMKIN species production rate function interface (e.g., QSS approximations)

Note that the CHEMKIN-II specification assumes a mixture of ideal gases.  If the `Pele` codes are built with `Eos_dir = Fuego`, the make system variable `Chemistry_Model` must be set to one of the models (subfolders) that exist in the `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models` folder. The repository currently provides (multiple) models for each of the following mixtures:
* air
* hydrogen
* methane
* ethane
* *n*-dodecane
* di-methyl ether
* ... more can be added (see below)


### NOTE: Non-ideal equations of state

`PelePhysics` currently supports a cubic EOS model: `Soave-Redlich-Kwong`.  It is built on top of the ideal gas models, and is selected by specifying its name as the `Eos_dir` during the build (the make system requires that both `Eos_dir` and `Chemistry_Model` be specified).  Any additional parameters required for the EOS (e.g., attractions, repulsions, critical states) are either included in the underlying FUEGO database used to generate the source file model implementation, or else are inferred from the input model data.

## Model generation procedures
This repository provides the tools necessary for generating new Pele-compatible combustion mechanisms.  To do this:
1. Ensure that the environment variable `PELE_PHYSICS_HOME` is set to the root folder of your local `PelePhysics` clone
2. Set up the required FUEGO variables by doing
   ```
   . ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/setup.sh` (for bash), or
   source ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/setup.csh` (for C-shell)
   ```
3. Check that all is setup correctly by reproducing one of the existing C++ model files.  For example, to generate the `LiDryer` model for H2 combustion:
```
     cd ${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/LiDryer  
     ./make-LiDryer.sh
```
   That script generates a C++ source file.  If this file "almost" matches the one in git, you're all set... By "almost" we mean minus reordering in the egtransetEPS, egtransetSIG, egtransetDIP, egtransetPOL, egtransetZROT and egtransetNLIN functions. If not, and there are no error messages, it might be that the git version needs to be updated.  Contact `MSDay@lbl.gov` to verify the correct actions to take.
4. To build and use a new CHEMKIN-II based model, make a new model folder in `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models`,  say `XXX`, and copy your CHEMKIN input files there.  Additionally copy `make-LiDryer.sh` from the `LiDryer` model folder there as a template, rename it to `make-XXX.sh`, and edit the filenames at the top.
5. The first step is to modify the chemistry input file to include a transport section. To do so:
```
${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models script_trans.py chem.inp trans.dat
```
 This should generate a text file, `TRANFILE_APPEND.txt`. Insert the entire contents of this file into the chemistry input file, just above REACTION section. Next, run your newly created script file, `make-XXX.sh`.  If all goes well, FUEGO will be called to generate a `chem.cpp` file that gets concatenated with a few others, into a result file, `XXX.cpp`, and cleans up its mess.
6. Add a `Make.package` text file in that model folder that will be included by the AMReX make system.  In most cases, this file will contain a single line, `CEXE_sources+=YYY.cpp` (see the other models for examples if there are auxiliary files in languages other than C++ to include in the build of this model).  
7. Finally, edit the `GNUmakefile` where you want to use this (in, e.g., `PeleC/Exec`) so that `Chemistry_Model` is `XXX`.  In `PeleC/Exec/Make.PeleC`, the model is expected to be in the folder `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/$(Chemistry_Model)`, and it is expected that the folder contains a `Make.package` file to include, so make sure things are where they need to be.
