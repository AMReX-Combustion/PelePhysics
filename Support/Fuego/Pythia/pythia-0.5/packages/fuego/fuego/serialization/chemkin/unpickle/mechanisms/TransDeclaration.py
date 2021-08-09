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
from .Declaration import Declaration


class TransDeclaration(Declaration):


    def __init__(self, symbol, locator=None):
        Declaration.__init__(self, locator)

        self.symbol = symbol
        self.phase = None
        self.composition = []

        return


    def __str__(self):
        str = "symbol='%s'" % self.symbol

        if self.phase:
            str += ", phase=%s" % self.phase

        if self.composition:
            str += ", composition=%s" % repr(self.composition)

        str += ", " + Declaration.__str__(self)

        return str


# version
__id__ = "$Id$"

# End of file
