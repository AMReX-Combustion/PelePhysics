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

from builtins import map
from builtins import object
from operator import mul
from functools import reduce


class PhaseSpace(object):


    def __init__(self, species, coefficients):
        self.species = species
        self.coefficients = coefficients
        return


    def __call__(self):
        concentrations = [ species.concentration() for species in self.species ]
        return reduce(mul, list(map(pow, concentrations, self.coefficients)))


# version
__id__ = "$Id$"

#
# End of file
