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
PACKAGE = serialization/mechanisms

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    Element.py \
    ElementSet.py \
    Entity.py \
    EntitySet.py \
    Mechanism.py \
    MechanismExceptions.py \
    NASA.py \
    Reaction.py \
    ReactionSet.py \
    Species.py \
    SpeciesSet.py \
    Thermodynamics.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
