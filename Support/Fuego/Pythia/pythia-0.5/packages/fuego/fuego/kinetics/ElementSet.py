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
    def element(self, symbol):
        return self.find(symbol)

    def __init__(self, elements):
        EntitySet.__init__(self)

        from .Element import Element

        for proxy in elements:
            element = Element(proxy.symbol, proxy.weight)
            self.insert(proxy.symbol, element)

        return


# version
__id__ = "$Id$"

# End of file
