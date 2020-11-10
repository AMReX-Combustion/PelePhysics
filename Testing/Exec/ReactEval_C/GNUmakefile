PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

TINY_PROFILE = TRUE

USE_CUDA   = FALSE

ifeq ($(USE_CUDA), TRUE)
    USE_OMP=FALSE
endif
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
      USE_KLU_PP = TRUE
  endif
<<<<<<< HEAD
  ifeq ($(USE_CUDA), TRUE)
  # Activates GPU version
  USE_CUDA_SUNDIALS_PP = FALSE
  endif
=======
>>>>>>> 9e95165c70bccecf7f25e11f610c6ee52724a013
else
  ifeq ($(USE_RK64_PP), TRUE)
  # Activates use of Hari explicit RK
  # will only work if USE_SUNDIALS_PP = FALSE
<<<<<<< HEAD
  TRANSPORT_CPP = TRUE
  endif
=======
  USE_RK64_PP = TRUE
>>>>>>> 9e95165c70bccecf7f25e11f610c6ee52724a013
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

LIBRARIES += -lnvToolsExt

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
