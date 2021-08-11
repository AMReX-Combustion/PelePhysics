# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = fuego
PACKAGE = serialization/f

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    FPickler.py \
    __init__.py


export:: export-package-python-modules

# version
# $Id$

# End of file
