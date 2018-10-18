# PelePhysics
*A respository of physics databases and implementation code for use with the `Pele` suite of of codes*

## FUEGO

This `PelePhysics` repository contains the FUEGO source code generation tool in order to support the inclusion
of chemical models specified in the CHEMKIN-II format.  Typically, such models include ASCII files specifying  
* chemistry (elements, chemical species, Arrhenius reaction parameters) in CHEMKIN-II format
* thermodynamics (NASA polynomial expressions)
* transport (input data for the CHEMKIN function `tranfit` to compute transport coefficients)
* ...optionally, there may be multiple source files to evaluate the species production rates using various reduced/QSS approximations

Note that the CHEMKIN-II specification assumes a mixture of ideal gases.  If the `Pele` codes are built with `Eos_dir = Fuego`, the variable `Chemistry_Model` must be set to one of the models existing in the `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models` folder (named after the folder containing the model files there).

The repository provides models for inert air, and reacting hydrogen, methane, ethane, dodecane and di-methyl ether, and more can be added, using the procedure below.  Note that although all provided models can be regenerated following these procedures, we include the resulting source files for convenience.  As a result, if the generation procedure is modified, the complete set of models will need to be updated manually.  There is currently no procure in place to keep the models up to date with FUEGO.

## Non-ideal equations of state

`PelePhysics` supports a number of EOS models, including the CHEMKIN-II implementation of ideal gas mixtures, a simple `GammaLaw` model, and a cubic model: `Soave-Redlich-Kwong`.  The cubic model is built on top of the ideal gas models, and is selected by specifying its name as the `Eos_dir` during the build (the `Chemistry_Model` must also be specified).  Any additional parameters (e.g., attractions, repulsions, critical states) are either included in the underlying FUEGO database used to generate the source file model implementation, or else are inferred from the input model data.

### Model generation procedures
1. Ensure that the environment variable `PELE_PHYSICS_HOME` is set to the root folder of your local `PelePhysics` clone containing this file.
2. Set up some fuego variables by doing `. ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/setup.sh` (for bash, or `source ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/setup.csh`)
3. Check that all is setup correctly by reproducing one of the existing .c files.  Say the LiDryer one...
```
     cd ${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/LiDryer  
     ./make-LiDryer.sh
```
   That script uses the CHEMKIN-II thermodynamics database files expected by the usual CHEMKIN programs as well as a chemistry input file modified to include a transport section (see below) and generates a c source file.  If this file "almost" matches the one in git, you're all set... By "almost" we mean minus reordering in the egtransetEPS, egtransetSIG, egtransetDIP, egtransetPOL, egtransetZROT and egtransetNLIN functions. If not, and there are no error messages, it might be that the git version needs to be updated.  Contact `MSDay@lbl.gov` to verify the correct actions to take.
5. To build and use a new CHEMKIN-II based model, make a new model folder in `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models`,  say `XXX`, and copy your CHEMKIN input files there.  Additionally copy `make-LiDryer.sh` from the `LiDryer` model folder there as a template, rename it to `make-XXX.sh`, and edit the filenames at the top.  You must start by modifying the chemistry input file to include a transport section. To do so, call the script located under `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models` labelled `script_trans.py` with your chem.inp and trans.dat as arguments (in that order). The script should return a text file labelled `TRANFILE_APPEND.txt`. Open this file, copy all the lines and paste them in your chemistry input file, say right above the REACTION section. Next, run the `make-XXX.sh` script file, and if all goes well, the software will generate a `chem.c` file that gets concatenated with a few others, puts everything into a result file, `YYY.c`, and cleans up its mess.
6. Add a `Make.package` text file in that model folder that will be included by the AMReX make system.  In most case, this file will contain a single line, `cEXE_sources+=YYY.c` (see the other models for examples if there are auxiliary files in languages other than c to include in the build of this model).  
7. Finally, edit the `GNUmakefile` where you want to use this (in, e.g., `PeleC/Exec`) so that `Chemistry_Model` is `XXX`.  In `PeleC/Exec/Make.PeleC`, the model is expected to be in the folder `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/$(Chemistry_Model)`, and it is expected that the folder contains a `Make.package` file to include, so make sure things are where they need to be.


Please see git@code.ornl.gov:Pele/PeleC/README.md for more information.

