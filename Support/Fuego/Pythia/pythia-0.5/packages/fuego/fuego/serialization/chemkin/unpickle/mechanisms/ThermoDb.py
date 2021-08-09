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


class ThermoDb(object):
    def all(self, flag=None):

        if flag is None:
            return self._all

        self._all = flag
        return

    def range(self, range=None):
        if range is None:
            return self._range

        self._range = range
        return

    def species(self, species):
        self._species.append(species)
        self._index[species.symbol] = species
        return

    def size(self):
        return len(self._species)

    def find(self, symbol=None):
        if symbol:
            return self._index.get(symbol)

        return self._species

    def __init__(self):
        self._all = 0
        self._range = ()

        self._index = {}
        self._species = []
        return


# version
__id__ = "$Id$"

# End of file
