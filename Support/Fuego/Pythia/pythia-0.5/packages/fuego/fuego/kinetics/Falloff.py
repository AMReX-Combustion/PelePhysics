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

from .ThirdBody import ThirdBody


class Falloff(ThirdBody):
    def update(self, T=None):

        # JES said: use the enhancement for the rate only
        # do not "double count" by using it as an overall factor as well

        enhancement = self.enhancementCalculator()
        self.phaseSpace = self.phaseSpaceCalculator()
        self.rate = self.rateCalculator(T, enhancement)
        return

    def __call__(self):
        return self.rate * self.phaseSpace


# version
__id__ = "$Id$"

#
# End of file
