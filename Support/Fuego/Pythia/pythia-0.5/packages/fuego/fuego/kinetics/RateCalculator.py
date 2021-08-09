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

from builtins import object

import pyre
from pyre.handbook.units.SI import kelvin


class RateCalculator(object):
    def update(self, T):
        raise NotImplementedError(
            "class '%s' should override 'update'" % (self.__class__.__name__)
        )

    def __init__(self, phaseSpaceCalculator, rateCalculator):
        self.phaseSpaceCalculator = phaseSpaceCalculator
        self.rateCalculator = rateCalculator

        self.rate = None
        self.enhancement = 1.0
        self.phaseSpace = None
        return

    def __call__(self):
        raise NotImplementedError(
            "class '%s' should override '__call__'" % (self.__class__.__name__)
        )

    def __str__(self):
        return self.__class__.__name__


# version
__id__ = "$Id$"

#
# End of file
