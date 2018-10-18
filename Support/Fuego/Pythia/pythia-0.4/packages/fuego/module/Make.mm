# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = fuego
PACKAGE = fuegomodule

include std-pythonmodule.def
include local.def

PROJ_CXX_SRCLIB = $(BLD_LIBDIR)/lib$(PROJECT).$(EXT_AR)

PROJ_SRCS = \
    bindings.cc \
    exceptions.cc \
    misc.cc


# version
# $Id$

# End of file
