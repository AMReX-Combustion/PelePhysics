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

from Arrhenius import Arrhenius


class ThirdBody(Arrhenius):


    def update(self, T=None):
        self.enhancement = self.enhancementCalculator()
        self.phaseSpace = self.phaseSpaceCalculator()
        self.rate = self.rateCalculator(T)
        return


    def __init__(self, phaseSpaceCalculator, enhancementCalculator, rateCalculator):
        Arrhenius.__init__(self, phaseSpaceCalculator, rateCalculator)
        self.enhancementCalculator = enhancementCalculator
        return


    def __call__(self):
        return self.enhancement * self.rate * self.phaseSpace


# version
__id__ = "$Id$"

#
# End of file
