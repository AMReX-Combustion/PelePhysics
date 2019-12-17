PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

USE_CUDA   = TRUE

TINY_PROFILE = TRUE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    := ../../..

#######################
# this flag activates the subcycling mode in the D/Cvode routines
DEFINES  += -DMOD_REACTOR

# use CVODE
#CVODE_LIB_DIR=$(PELE_PHYSICS_HOME)/ThirdParty/sundials/instdir/lib/
USE_SUNDIALS_PP = TRUE

USE_CUDA_CVODE = TRUE

#######################
ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  Chemistry_Model = drm19
  #Chemistry_Model = dodecane_wang
  #Chemistry_Model = LiDryer
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
