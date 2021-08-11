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

PROJECT = weaver
PACKAGE = mills

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    BlockComments.py \
    BlockMill.py \
    CMill.py \
    CommentingStrategy.py \
    CshMill.py \
    CxxMill.py \
    Fortran77Mill.py \
    Fortran90Mill.py \
    HTMLMill.py \
    Indenter.py \
    LineComments.py \
    LineMill.py \
    MakeMill.py \
    Mill.py \
    MillRegistrar.py \
    PerlMill.py \
    PythonMill.py \
    ShMill.py \
    Stationery.py \
    TeXMill.py \
    XMLMill.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
