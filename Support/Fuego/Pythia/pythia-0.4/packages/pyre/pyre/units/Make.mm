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
PACKAGE = units


#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    SI.py \
    area.py \
    density.py \
    energy.py \
    force.py \
    length.py \
    mass.py \
    power.py \
    pressure.py \
    speed.py \
    substance.py \
    temperature.py \
    time.py \
    unit.py \
    unitparser.py \
    volume.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
