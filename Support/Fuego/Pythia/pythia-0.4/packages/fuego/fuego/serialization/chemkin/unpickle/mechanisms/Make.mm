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
PACKAGE = serialization/chemkin/unpickle/mechanisms

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    Declaration.py \
    ElementDb.py \
    ElementDeclaration.py \
    ExternalThermo.py \
    Mechanism.py \
    MechanismExceptions.py \
    ReactionDb.py \
    ReactionDeclaration.py \
    SpeciesDb.py \
    SpeciesDeclaration.py \
    ThermoDb.py \
    ThermoDeclaration.py \
    TransDb.py \
    TransDeclaration.py\
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
