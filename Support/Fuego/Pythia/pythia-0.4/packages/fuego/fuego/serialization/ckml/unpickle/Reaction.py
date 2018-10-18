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

import pyre.util.bool
from AbstractNode import AbstractNode
 

class Reaction(AbstractNode):


    tag = "reaction"


    def notify(self, parent):
        parent.onReaction(self._reaction)
        return


    def onDuplicate(self):
        self._reaction.duplicate = 1
        return


    def onReactionUnits(self, units):
        self._reaction.units = units
        return


    def onReagents(self, reactants, products):
        self._reaction.reactants = reactants
        self._reaction.products = products
        return


    def onReactionRate(self, rate):
        self._reaction.efficiencies = rate.efficiencies

        self._reaction.arrhenius = rate.arrhenius
        self._reaction.rev = rate.rev
        self._reaction.lt = rate.lt
        self._reaction.rlt = rate.rlt

        self._reaction.low = rate.low
        self._reaction.sri = rate.sri
        self._reaction.troe = rate.troe

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        symbol = int(attributes["id"])

        falloff = pyre.util.bool.bool(attributes.get("falloff", "false"))
        thirdBody = ()
        species = attributes.get("thirdBody")
        if species:
            if species == "mixture":
                thirdBody = ("<mixture>", 1)
            else:
                thirdBody = (species, 1)
            
        reversible = pyre.util.bool.bool(attributes.get("reversible", "false"))

        locator = self.documentNode().locator()

        reaction = self.documentNode().mechanism().newReaction(symbol, locator)
        reaction.falloff = falloff
        reaction.thirdBody = thirdBody
        reaction.reversible = reversible

        self._reaction = reaction

        return
            

# version
__id__ = "$Id$"

#  End of file 
