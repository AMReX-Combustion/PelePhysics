#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

from __future__ import absolute_import
from .AbstractNode import AbstractNode
 

class Species(AbstractNode):


    tag = "species"


    def notify(self, parent):
        parent.onSpecies(self._species)
        return


    def onComposition(self, composition):
        self._species.composition = composition
        return


    def onThermo(self, thermo):
        for p in thermo:
            type, lowT, highT, parameters = p
            self._species.thermalParametrization(
                type, lowT, highT, self.documentNode().locator(), parameters)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        symbol = attributes["id"]
        phase = attributes.get("phase")
        locator = self.documentNode().locator()

        species = self.documentNode().mechanism().newSpecies(symbol, locator)
        species.phase = phase
        self._species = species
        return
            

# version
__id__ = "$Id$"

#  End of file 
