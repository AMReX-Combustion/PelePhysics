PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_MPI    = TRUE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    := ../../..

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  #Chemistry_Model = chem-H
  Chemistry_Model = LiDryer
  #Chemistry_Model = BurkeDryer_mod
  #Chemistry_Model = BurkeDryer
  #Chemistry_Model = air
  #Chemistry_Model = grimech12
  #Chemistry_Model = grimech30
  #Chemistry_Model = grimech30-noArN
  #Chemistry_Model = drm19
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
