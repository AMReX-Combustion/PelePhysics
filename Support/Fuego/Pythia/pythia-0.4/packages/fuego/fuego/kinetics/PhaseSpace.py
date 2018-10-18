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

from operator import mul


class PhaseSpace:


    def __init__(self, species, coefficients):
        self.species = species
        self.coefficients = coefficients
        return


    def __call__(self):
        concentrations = [ species.concentration() for species in self.species ]
        return reduce(mul, map(pow, concentrations, self.coefficients))


# version
__id__ = "$Id$"

#
# End of file
