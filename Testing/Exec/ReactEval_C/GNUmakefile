PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

#TINY_PROFILE = TRUE

USE_CUDA   = TRUE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    := ../../..

#######################
# this flag activates the subcycling mode in the D/Cvode routines
DEFINES  += -DMOD_REACTOR

#######################
# ODE solver OPTIONS: DVODE (default) / SUNDIALS / RK explicit
#######################
# Activates use of SUNDIALS: CVODE (default) / ARKODE
USE_SUNDIALS_PP = TRUE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  # provide location of sundials lib if needed
  #SUNDIALS_LIB_DIR=$(PELE_PHYSICS_HOME)/ThirdParty/sundials/instdir/lib/
  # use ARKODE: if FALSE then CVODE is used
  USE_ARKODE_PP = FALSE
  # use KLU sparse features -- only useful if CVODE is used on CPU
  USE_KLU_PP = FALSE
  ifeq ($(USE_KLU_PP), TRUE)
    # provide location of KLU lib if needed
    #SUITESPARSE_DIR=$(PELE_PHYSICS_HOME)/ThirdParty/SuiteSparse/
  endif
  # Activates GPU version
  USE_CUDA_SUNDIALS_PP = TRUE
else
  # Activates use of Hari explicit RK
  # will only work if USE_SUNDIALS_PP = FALSE
  USE_RK64_PP = FALSE
endif

#######################
ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  Chemistry_Model = drm19
  Reactions_dir   = Fuego
  Transport_dir   = Simple
else
  Eos_dir       = GammaLaw
  Reactions_dir = Null
  Transport_dir = Constant
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
