PELE_PHYSICS_HOME    := ../../..


PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_CUDA   = FALSE
USE_CUDA   = TRUE

TINY_PROFILE = FALSE

USE_SUNDIALS_PP = TRUE
ifeq ($(USE_CUDA), TRUE)
  USE_KLU_PP = FALSE
else
  USE_F90_PP = FALSE
  USE_KLU_PP = FALSE
endif
USE_KLU = $(USE_KLU_PP) # sort out keywords later...

Eos_Model       = Fuego
#Chemistry_Model = LiDryer
Chemistry_Model = drm19
#Chemistry_Model = grimech12
Transport_Model = Simple

CEXE_sources += main.cpp
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
