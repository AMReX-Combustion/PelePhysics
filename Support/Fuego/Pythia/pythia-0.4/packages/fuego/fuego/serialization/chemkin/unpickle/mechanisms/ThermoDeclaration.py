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

from Declaration import Declaration


class ThermoDeclaration(Declaration):


    def __init__(self, symbol, locator=None):
        Declaration.__init__(self, locator)

        self.symbol = symbol
        self.phase = None
        self.composition = []

        self.range = None
        self.lowParameters = []
        self.highParameters = []

        return


    def __str__(self):
        str = "symbol='%s'" % self.symbol

        if self.phase:
            str += ", phase=%s" % self.phase

        if self.composition:
            str += ", composition=%s" % `self.composition`

        if self.range:
            str += ", parametrization range=%s" % `self.range`

        if self.lowParameters:
            str += ", low=%s" % `self.lowParameters`

        if self.highParameters:
            str += ", high=%s" % `self.highParameters`

        str += ", " + Declaration.__str__(self)

        return str


# version
__id__ = "$Id$"

# End of file
