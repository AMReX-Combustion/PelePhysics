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
class QssSpeciesDb(object):


    def qss_species(self, qss_species):
        self._qss_species.append(qss_species)
        self._qss_index[qss_species.symbol] = qss_species
        return


    def size(self):
        return len(self._qss_species)


    def find(self, symbol=None):
        if symbol:
            return self._qss_index.get(symbol)

        return self._qss_species


    def __init__(self):
        self._qss_index = {}
        self._qss_species = []
        return


# version
__id__ = "$Id$"

# End of file
