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

PROJECT = journal
PACKAGE = _journalmodule
MODULE = _journal

include std-pythonmodule.def

PROJ_CXX_SRCLIB = -ljournal

PROJ_SRCS = \
    ProxyCategory.cc \
    ProxyDevice.cc \
    ProxyIndex.cc \
    bindings.cc \
    exceptions.cc \
    journal.cc \
    misc.cc \


# version
# $Id$

# End of file
