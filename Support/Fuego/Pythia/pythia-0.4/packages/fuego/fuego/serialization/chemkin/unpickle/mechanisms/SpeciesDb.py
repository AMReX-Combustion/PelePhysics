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


class SpeciesDb:


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
        self._index = {}
        self._species = []
        return


# version
__id__ = "$Id$"

# End of file
