# AMReX
DIM        = 2
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = FALSE
DEBUG      = TRUE

# Compiler
COMP	     = gnu
FCOMP	     = gfortran
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE

# PelePhysics
FUEGO_GAS  = TRUE
USE_SUNDIALS_PP = FALSE

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  Chemistry_Model = air
  Reactions_dir   = Null
  Transport_dir   = Simple
else
  Eos_dir         = GammaLaw
  Chemistry_Model = Null
  Reactions_dir   = Null
  Transport_dir   = Constant
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
