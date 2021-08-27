# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = FALSE
DEBUG      = FALSE

# Compiler
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE

Eos_Model       = GammaLaw
Chemistry_Model = Null
Reactions_dir   = null
Transport_Model = Constant

CEXE_sources += main.cpp
Bpack   := $(PELE_PHYSICS_HOME)/Utility/Make.package
Blocs   := . $(PELE_PHYSICS_HOME)/Utility

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
