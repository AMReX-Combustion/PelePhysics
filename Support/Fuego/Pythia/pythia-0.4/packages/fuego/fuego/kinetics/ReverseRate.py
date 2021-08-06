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
from __future__ import division
import pyre
from pyre.handbook.units.SI import kelvin


from .RateCalculator import RateCalculator


class ReverseRate(RateCalculator):


    def update(self, T=None):

        self.phaseSpace = self.phaseSpaceCalculator()
        self.enhancement = self.rateCalculator.enhancement
        self.rate = self.rateCalculator.rate / self.equilibriumCalculator(T)

        return


    def temperature(self, T):
        return


    def __init__(self, phaseSpaceCalculator, forwardRate, equilibriumCalculator):
        RateCalculator.__init__(self, phaseSpaceCalculator, forwardRate)

        self.equilibriumCalculator = equilibriumCalculator

        return


    def __call__(self):
        return self.enhancement * self.rate * self.phaseSpace


# version
__id__ = "$Id$"

#
# End of file
