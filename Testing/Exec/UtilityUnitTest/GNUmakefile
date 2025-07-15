# AMReX
DIM          = 2
PRECISION    = DOUBLE
PROFILE      = FALSE
TINY_PROFILE = FALSE
VERBOSE      = FALSE
DEBUG        = FALSE

# Compiler
COMP       = gnu
FCOMP      = gfortran
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE

PELE_PHYSICS_HOME ?= ../../..
AMREX_HOME        ?= $(abspath $(lastword $(PELE_PHYSICS_HOME)/Submodules/amrex))
UTILITY_HOME      ?= $(PELE_PHYSICS_HOME)/Source/Utility

EBASE = PelePhysicsUnitTest
BL_NO_FORT = TRUE

##

Bpack   := ./Make.package
Blocs   := .

PP_SRC_HOME = $(PELE_PHYSICS_HOME)/Source
EXTERN_CORE       += $(PP_SRC_HOME)
INCLUDE_LOCATIONS += $(PP_SRC_HOME)
VPATH_LOCATIONS   += $(PP_SRC_HOME)
Bpack             += $(PP_SRC_HOME)/Make.package
Blocs             += $(PP_SRC_HOME)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

INCLUDE_LOCATIONS += $(UTILITY_HOME)
VPATH_LOCATIONS   += $(UTILITY_HOME)
#include $(UTILITY_HOME)/Diagnostics/Make.package
#include $(UTILITY_HOME)/Filter/Make.package
#include $(UTILITY_HOME)/PMF/Make.package # Needs mechanism.H but shouldn't
#include $(UTILITY_HOME)/PltFileManager/Make.package
#include $(UTILITY_HOME)/TurbInflow/Make.package
#include $(UTILITY_HOME)/BlackBoxFunction/Make.package
include $(UTILITY_HOME)/Utilities/Make.package

Pdirs   := Base AmrCore Boundary
Bpack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))
include $(Bpack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

PP_SRC_HOME = $(PELE_PHYSICS_HOME)/Source
EXTERN_CORE       += $(PP_SRC_HOME)
INCLUDE_LOCATIONS += $(PP_SRC_HOME)
VPATH_LOCATIONS   += $(PP_SRC_HOME)
Bpack             += $(PP_SRC_HOME)/Make.package
Blocs             += $(PP_SRC_HOME)

# job_info support
CEXE_sources += AMReX_buildInfo.cpp
CEXE_headers += $(AMREX_HOME)/Tools/C_scripts/buildInfo.H
INCLUDE_LOCATIONS +=  $(AMREX_HOME)/Tools/C_scripts

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) AMReX_buildInfo.cpp

AMReX_buildInfo.cpp:
	$(AMREX_HOME)/Tools/C_scripts/makebuildinfo_C.py \
          --amrex_home "$(AMREX_HOME)" \
          --COMP "$(COMP)" --COMP_VERSION "$(COMP_VERSION)" \
          --FCOMP "$(FCOMP)" --FCOMP_VERSION "$(FCOMP_VERSION)" \
          --MODULES "$(MNAMES)" \
          --GIT "$(PELE_PHYSICS_HOME) $(AMREX_HOME) $(SUNDIALS_HOME)"
