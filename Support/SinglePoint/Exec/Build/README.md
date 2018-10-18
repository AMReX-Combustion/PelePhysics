# SinglePoint

Evaluation of transport coefficients at a single point, using the "Constant" or "Simple" model.  The
"Simple" model is a stripped down version of the mixtured-average transport model evaluated by 
EGLib (see http://www.cmap.polytechnique.fr/www.eglib), but is suitable for public distribution.
In the future, a procedure will be available to link this library directly with EGLib routines.

## Build/running the demo

This folder provides a capability for building a standalone library for evaluating the transport
models, and, optionally, a driver program to demonstrate how to call the function for a single data point.

The GNUmakefile in this folder uses the AMReX makefile system to build a shared library, libSP.XXX.so, 
(XXX is a build-specific string encoding compiler and physics options choices, allowing multiple versions
to co-exist).  Setting the makefile variable, BUILD_EXEC=TRUE, will link the library with `main.cpp`
to demonstrate usage.  The choices of EOS and transport models are set with makefile variables in the 
GNUmakefile.

Currently, this standalone setup supports "gamma law" gases with constant transport coefficients,
or "fuego" gases (which are mixtures of ideal gases, specified originally in CHEMKIN-II format, and
converted to compiled form by the LBNL version of the FUEGO app), with "Simple" transport.  The
FUEGO app, included in this distribution, is discussed in the top-level README file.

If BUILD_EXEC=TRUE, and executable named SPXXX.ex is built. Note that both transport model implementations 
driven by this demo expect a fortran namelist file in the folder where the executable is run.  The namelist file
may be empty, and its name is passed to the fort_data_init prior to evaluating any properties. The contents
of the namelist file depend on the model, and samples are provided for the "Constant" and "Simple" models.

* The fuego-based model requires that a "Chemistry_Model" be specified in the GNUmakefile.  It does not
currently make use of any namelisted variables.

* The constant transport coefficient model looks for a namelist called `extern`, which may contain
values for `const_conductivity`, `const_viscosity`, `const_bulk_viscosity`, `const_diffusivity`.  Default
values for these variables are automatically set in a generated source file, extern.f90, during the
build.

You MUST have a namelist file present to run this demo driver, so you might copy/link one of the
provided samples to "probin" (or set the appropriate name of your file in main.cpp or your own calling routine).

See the GNUmakefile for additional build details, examples.

If you run into difficulties, contact Marc Day at MSDay@lbl.gov, or file a GitLab "issue" for support.
