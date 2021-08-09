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


class ElementDeclaration(Declaration):
    def __init__(self, symbol, weight=None, locator=None):
        Declaration.__init__(self, locator)

        self.symbol = symbol
        self.weight = weight

        return

    def __str__(self):
        str = "symbol='%s'" % self.symbol
        if self.weight:
            str += "weight=%g" % self.weight

        str += ", " + Declaration.__str__(self)

        return str


# version
__id__ = "$Id$"

# End of file
