PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
#FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

FUEGO_GAS  = TRUE

USE_CUDA   = TRUE

TINY_PROFILE = TRUE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    := ../../..

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  #Chemistry_Model = LiDryer
  Chemistry_Model = drm19
  #Chemistry_Model = grimech12
  #Chemistry_Model = grimech30
  #Chemistry_Model = dodecane_wang
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
