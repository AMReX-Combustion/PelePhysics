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

from .EntitySet import EntitySet


class ElementSet(EntitySet):
    def element(self, symbol, atomicWeight=None, locator=None):
        from .Element import Element

        element = Element(self.size(), symbol, atomicWeight, locator)
        self.insert(symbol, element)
        return element


# version
__id__ = "$Id$"

# End of file
