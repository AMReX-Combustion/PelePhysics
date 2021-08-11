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

from builtins import object


class ExternalThermo(object):
    def element(self, symbol=None):
        return symbol

    def species(self, symbol=None):
        return symbol

    # thermo

    def declareThermalProperties(self, species):
        symbol = species.symbol
        duplicate = self._thermoDb.find(symbol)

        self._thermoDb.species(species)

        if duplicate:
            raise self.DuplicateThermo(symbol)

        return

    def thermalProperties(self, species=None):
        prop = self._thermoDb.find(species)
        if not prop:
            return self._externalDb.find(species)

    def thermalDatabase(self):
        return self._thermoDb

    def thermoAll(self):
        return self._thermoDb.all(1)

    def thermoDone(self):
        return

    def thermoRange(self, range=None):
        return self._thermoDb.range(range)

    # other methods

    def __init__(self, source):

        self._source = source

        from .ThermoDb import ThermoDb

        self._thermoDb = ThermoDb()

        from .ThermoDeclaration import ThermoDeclaration

        self.thermoFactory = ThermoDeclaration

        return


# version
__id__ = "$Id$"

# End of file
