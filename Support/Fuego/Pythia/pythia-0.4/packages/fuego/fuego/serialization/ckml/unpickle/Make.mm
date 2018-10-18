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
PACKAGE = serialization/ckml/unpickle

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    AbstractNode.py \
    Arrhenius.py \
    Atom.py \
    CkmlParser.py \
    Composition.py \
    Document.py \
    Duplicate.py \
    Efficiencies.py \
    Element.py \
    Enhancement.py \
    Kinetics.py \
    LandauTeller.py \
    LowPressure.py \
    NASA.py \
    NASA_a1.py \
    NASA_a2.py \
    NASA_a3.py \
    NASA_a4.py \
    NASA_a5.py \
    NASA_a6.py \
    NASA_a7.py \
    Product.py \
    Reactant.py \
    Reaction.py \
    ReactionRate.py \
    ReactionUnits.py \
    Reagents.py \
    Reverse.py \
    ReverseLandauTeller.py \
    SRI.py \
    Species.py \
    TROE.py \
    Thermo.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
