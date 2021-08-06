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
import pyre
from pyre.handbook.units.SI import kelvin


from .RateCalculator import RateCalculator


class Arrhenius(RateCalculator):


    def update(self, T):
        self.phaseSpace = self.phaseSpaceCalculator()
        self.rate = self.rateCalculator(T)
        return


    def __call__(self):
        return self.rate * self.phaseSpace


# version
__id__ = "$Id$"

#
# End of file
