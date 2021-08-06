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
from .Entity import Entity


class Element(Entity):


    def __init__(self, id, symbol, weight=None, locator=None):
        Entity.__init__(self, id, locator)
        
        self.symbol = symbol
        self.weight = weight

        return


    def __str__(self):
        str = "symbol='%s'" % self.symbol
        if self.weight:
            str += "weight=%g" % self.weight

        str += ", source=" + Entity.__str__(self)

        return str


    __slots__ = ("symbol", "weight")


# version
__id__ = "$Id$"

# End of file
