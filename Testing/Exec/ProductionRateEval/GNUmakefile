PELE_PHYSICS_HOME    := ../../..


PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_CUDA   = TRUE
USE_CUDA   = FALSE

TINY_PROFILE = FALSE

#######################
USE_SUNDIALS_PP = TRUE
USE_SUNDIALS_PP = FALSE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  # use KLU sparse features -- only useful if CVODE is used on CPU
  USE_KLU_PP = FALSE
  # Activates GPU version
  USE_CUDA_SUNDIALS_PP = TRUE
endif
#######################

Eos_dir         = Fuego
Chemistry_Model = LiDryer
Reactions_dir   = Fuego
Transport_dir   = Simple

CEXE_sources += main.cpp
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
