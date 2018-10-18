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
PACKAGE = geometry/unpickle/pml


#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    AbstractNode.py \
    Angle.py \
    Binary.py \
    Block.py \
    Composition.py \
    Cone.py \
    Cylinder.py \
    Difference.py \
    Dilation.py \
    Document.py \
    GeneralizedCone.py \
    Geometry.py \
    Intersection.py \
    Prism.py \
    Pyramid.py \
    Reflection.py \
    Reversal.py \
    Rotation.py \
    Scale.py \
    Sphere.py \
    Torus.py \
    Transformation.py \
    Translation.py \
    Union.py \
    Vector.py \
    __init__.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
