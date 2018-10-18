# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def

PROJECT = journal
PACKAGE = libjournal

PROJ_DLL = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_SO)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)/$(PACKAGE)
PROJ_CLEAN += $(PROJ_INCDIR) $(PROJ_DLL)

PROJ_SRCS = \
    Category.cc \
    Console.cc \
    Debug.cc \
    DefaultRenderer.cc \
    Device.cc \
    Diagnostic.cc \
    Entry.cc \
    Error.cc \
    Firewall.cc \
    Index.cc \
    Info.cc \
    Journal.cc \
    LocalCategory.cc \
    LocalIndex.cc \
    Renderer.cc \
    StreamDevice.cc \
    Warning.cc \
    debuginfo.cc \
    firewall.cc

all: $(PROJ_DLL) export

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build the shared object

$(PROJ_DLL): product_dirs $(PROJ_OBJS)
	$(CXX) -o $(PROJ_DLL) $(PROJ_OBJS) $(LCXXFLAGS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
export:: export-headers export-libraries

EXPORT_HEADERS = \
    Category.h \
    Category.icc \
    Debug.h \
    Diagnostic.h \
    Diagnostic.icc \
    Error.h \
    Firewall.h \
    Info.h \
    Warning.h \
    journal.h \
    manipulators.h \
    manip-decls.h \
    manip-decls.icc \
    debuginfo.h \
    firewall.h \
    macros.h


EXPORT_LIBS = $(PROJ_DLL)


# version
# $Id$

#
# End of file
