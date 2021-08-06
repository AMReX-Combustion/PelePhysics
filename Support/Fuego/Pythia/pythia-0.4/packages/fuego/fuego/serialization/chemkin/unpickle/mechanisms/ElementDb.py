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
class ElementDb(object):


    def element(self, element):
        self._elements.append(element)
        self._index[element.symbol] = element
        return


    def size(self):
        return len(self._elements)


    def find(self, symbol=None):
        if symbol:
            return self._index.get(symbol)

        return self._elements


    def __init__(self):
        self._index = {}
        self._elements = []
        return


# version
__id__ = "$Id$"

# End of file
