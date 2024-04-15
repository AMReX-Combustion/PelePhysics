AMREX_HOME   ?= ../../Submodules/amrex

# Local
EBASE = PeleTurb

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Bpack   += ./Make.package

Pdirs   := Base Boundary
Bpack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Blocs   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))
include $(Bpack)
INCLUDE_LOCATIONS += $(Blocs)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
