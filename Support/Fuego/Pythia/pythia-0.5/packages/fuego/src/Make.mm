# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def

PROJECT = fuego
PACKAGE = libfuego

PROJ_CXX_LIB = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_LIB)
PROJ_DLL = $(PACKAGE).$(EXT_SO)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)/$(PACKAGE)
PROJ_LIBRARIES = 
PROJ_CLEAN += $(PROJ_INCDIR)

PROJ_SRCS = \
    hello.cc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build the library

all: $(PROJ_CXX_LIB)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build the shared object

$(PROJ_DLL): $(PROJ_OBJS)
	$(CXX) -o $(PROJ_DLL) $(PROJ_OBJS) $(LCXX_SOFLAGS) $(LCXXFLAGS) $(PROJ_LIBRARIES)
	$(CP_F) $(PROJ_DLL) $(BLD_LIBDIR)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# export

export:: export-headers export-binaries

EXPORT_HEADERS = \
    hello.h

EXPORT_BINS = $(PROJ_DLL)


# version
# $Id$

#
# End of file
