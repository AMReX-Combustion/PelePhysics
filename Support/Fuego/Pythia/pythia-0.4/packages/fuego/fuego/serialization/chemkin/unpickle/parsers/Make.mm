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
PACKAGE = serialization/chemkin/unpickle/parsers

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    BaseParser.py \
    ChemkinParser.py \
    Elements.py \
    Parameters.py \
    QssSpecies.py \
    Reactions.py \
    Species.py \
    Thermo.py \
    ThermoDatabaseParser.py \
    Trans.py \
    TransDatabaseParser.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
