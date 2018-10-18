# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = pyre
PACKAGE = properties


#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    Binary.py \
    Bool.py \
    Choice.py \
    Dimensional.py \
    Float.py \
    Greater.py \
    Integer.py \
    Less.py \
    List.py \
    Property.py \
    Range.py \
    Sequence.py \
    String.py \
    Unary.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
