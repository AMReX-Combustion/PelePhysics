#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from __future__ import absolute_import
def mixture(speciesSet):
    from .Mixture import Mixture
    return Mixture(speciesSet)


def mechanism(mix, proxy):
    from .Mechanism import Mechanism
    return Mechanism(mix, proxy)


def element(symbol, atomicWeight=None):
    from .Element import Element
    return Element(symbol, atomicWeight)


def species(symbol, phase, composition, thermalProperties):
    from .Species import Species
    return Species(symbol, phase, composition, thermalProperties)


# version
__id__ = "$Id$"

# End of file
