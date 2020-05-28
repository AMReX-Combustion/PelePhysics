PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = TRUE

FUEGO_GAS  = TRUE

#Determine the type of test ?
REGTEST_CP = FALSE
ifeq ($(REGTEST_CP), TRUE)
  DEFINES  += -DREGTEST_CP
endif
REGTEST_CV = FALSE
ifeq ($(REGTEST_CV), TRUE)
  DEFINES  += -DREGTEST_CV
endif
REGTEST_EY2T = FALSE
ifeq ($(REGTEST_EY2T), TRUE)
  DEFINES  += -DREGTEST_EY2T
endif

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    := ../../..

ifeq ($(FUEGO_GAS), TRUE)
  Eos_dir         = Fuego
  Chemistry_Model = LiDryer
  Reactions_dir   = Fuego
  Transport_dir   = Simple
else
  Eos_dir         = GammaLaw
  Chemistry_Model = Null
  Reactions_dir   = Fuego
  Transport_dir   = Constant
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
