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

PROJECT = fuego
PACKAGE = kinetics


#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    Arrhenius.py \
    EfficiencyEnhancement.py \
    Element.py \
    ElementSet.py \
    EntitySet.py \
    Equilibrium.py \
    Falloff.py \
    Mechanism.py \
    Mixture.py \
    PhaseSpace.py \
    RateCalculator.py \
    Reaction.py \
    ReactionSet.py \
    ReverseRate.py \
    Species.py \
    ThirdBody.py \
    ThirdBodyEnhancement.py \
    rateFactory.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
