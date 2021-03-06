# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = FALSE
DEBUG      = FALSE

# Compiler
COMP	     = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE

# PelePhysics
FUEGO_GAS  = TRUE
TINY_PROFILE = FALSE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= ../../..

#######################
# this flag activates the subcycling mode in the D/Cvode routines
DEFINES  += -DMOD_REACTOR

#######################
# ODE solver OPTIONS: DVODE (default) / SUNDIALS / RK explicit
#######################
# Activates use of SUNDIALS: CVODE (default) / ARKODE
USE_SUNDIALS_PP = TRUE
USE_RK64_PP = FALSE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  USE_ARKODE_PP = FALSE
  ifeq ($(USE_CUDA), TRUE)
      DEFINES += -DCVODE_BOXINTEG
      USE_KLU_PP = FALSE
  else
      ifeq ($(USE_HIP), TRUE)
          DEFINES += -DCVODE_BOXINTEG
          USE_KLU_PP = FALSE
      else
         USE_KLU_PP = TRUE
      endif
  endif
else
  USE_RK64_PP = TRUE
  ifeq ($(USE_RK64_PP), TRUE)
  # Activates use of Hari explicit RK
  # will only work if USE_SUNDIALS_PP = FALSE
  TRANSPORT_CPP = TRUE
  endif
endif

#######################
ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  Chemistry_Model = drm19
  Reactions_dir   = Fuego
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
