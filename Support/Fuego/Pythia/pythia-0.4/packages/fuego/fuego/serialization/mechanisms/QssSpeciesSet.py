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


class QssSpeciesSet(EntitySet):


    def qss_species(self, symbol, locator=None):
        from .QssSpecies import QssSpecies
        qss_species = QssSpecies(self.size(), symbol, locator)
        self.insert(symbol, qss_species)
        return qss_species



# version
__id__ = "$Id$"

# End of file
