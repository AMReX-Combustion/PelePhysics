PELE_PHYSICS_HOME    := ../../..


PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_CUDA   = TRUE
USE_CUDA   = FALSE

TINY_PROFILE = FALSE

ifeq ($(USE_CUDA), TRUE)
  PELE_USE_KLU = FALSE
else
  PELE_USE_KLU = FALSE
endif
USE_KLU = $(PELE_USE_KLU) # sort out keywords later...

Eos_Model       = Fuego
#Chemistry_Model = LiDryer
Chemistry_Model = drm19
#Chemistry_Model = grimech12
Transport_Model = Simple

CEXE_sources += main.cpp
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
