PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 1

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

BUILD_EXEC = TRUE
BUILD_EXEC = FALSE

# If FUEGO_GAS=TRUE, a Chemistry_Model must be defined, and Transport_Model may be Constant or Simple
# If FUEGO_GAS=FALSE, a GammaLaw gas is assumed, and the Transport_Model must be Constant
FUEGO_GAS = FALSE
FUEGO_GAS = TRUE

ifeq ($(FUEGO_GAS), TRUE)
  #Chemistry_Model = chem-H
  Chemistry_Model = LiDryer
  #Chemistry_Model = BurkeDryer_mod
  #Chemistry_Model = BurkeDryer
  #Chemistry_Model = air
  #Chemistry_Model = grimech12
  #Chemistry_Model = grimech30
  #Chemistry_Model = grimech30-noArN
  #Chemistry_Model = drm19
  Transport_dir   = Constant
  Transport_dir   = Simple
else
  Transport_dir = Constant
endif

ifeq ($(BUILD_EXEC), TRUE)
  CEXE_sources += main.cpp

  all: my_exec
endif

include ../Make.SinglePoint

my_exec: $(SPLIB) $(executable) 
	$(SILENT) $(RM) AMReX_buildInfo.cpp
	@echo SUCCESS
