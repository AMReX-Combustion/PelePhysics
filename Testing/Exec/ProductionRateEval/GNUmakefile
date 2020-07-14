PELE_PHYSICS_HOME    := ../../..


PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_CUDA   = FALSE
USE_CUDA   = TRUE

TINY_PROFILE = FALSE

TRANSPORT_CPP = TRUE

ifeq ($(USE_CUDA), TRUE)
  USE_SUNDIALS_PP = TRUE
  USE_KLU_PP = FALSE
  USE_CUDA_SUNDIALS_PP = TRUE
else
  USE_SUNDIALS_PP = FALSE
  USE_KLU_PP = FALSE
  USE_CUDA_SUNDIALS_PP = FALSE
endif
#######################

Eos_dir         = Fuego
Chemistry_Model = LiDryer
Reactions_dir   = Fuego
Transport_dir   = Simple

CEXE_sources += main.cpp
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
