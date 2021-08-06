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


class ReactionSet(EntitySet):


    def reaction(self, id, locator):
        from .Reaction import Reaction
        reaction = Reaction(id, locator)
        self.insert(reaction, reaction)
        return reaction


    def find(self, symbol=None, id=None):
        if id is not None:
            return self._entities[id]

        if not symbol:
            return self._entities

        candidates = []
        for reaction in self._entities:
            for species, coefficient in reaction.reactants:
                if species == symbol:
                    candidates.append(reaction)
                continue
            for species, coefficient in reaction.products:
                if species == symbol:
                    candidates.append(reaction)
                continue
                
        return candidates


# version
__id__ = "$Id$"

# End of file
