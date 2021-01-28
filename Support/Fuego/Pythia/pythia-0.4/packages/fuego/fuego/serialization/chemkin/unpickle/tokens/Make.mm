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
PACKAGE = serialization/chemkin/unpickle/tokens

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    ArrheniusCoefficients.py \
    Comments.py \
    ElementName.py \
    ElementSection.py \
    EndSection.py \
    ParameterList.py \
    Reaction.py \
    ReactionDuplicate.py \
    ReactionEfficiencies.py \
    ReactionFORD.py \
    ReactionHV.py \
    ReactionLOW.py \
    ReactionLT.py \
    ReactionRLT.py \
    ReactionReverse.py \
    ReactionSRI.py \
    ReactionSection.py \
    ReactionTROE.py \
    ReactionUnitsA.py \
    ReactionUnitsE.py \
    RegularExpressions.py \
    SpeciesName.py \
    QssSpeciesName.py \
    SpeciesSection.py \
    QssSpeciesSection.py \
    TemperatureRange.py \
    ThermoLine.py \
    ThermoSection.py \
    Token.py \
    Whitespace.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
