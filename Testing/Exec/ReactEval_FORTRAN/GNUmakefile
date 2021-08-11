PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

TINY_PROFILE = FALSE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    := ../../..

#Testing F90 routines
USE_F90_PP = TRUE

#######################
ifeq ($(FUEGO_GAS), TRUE)
  Eos_Model       = Fuego
  Chemistry_Model = LiDryer
  Reactions_dir   = Fuego
  Transport_Model = Simple
else
  Eos_Model       = GammaLaw
  Reactions_dir   = Null
  Transport_Model = Constant
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
