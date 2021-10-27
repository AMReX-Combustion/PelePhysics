# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= ../../..

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

# PelePhysics
TINY_PROFILE = FALSE

Eos_Model       = Fuego
Chemistry_Model = drm19
Transport_Model = Simple

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
