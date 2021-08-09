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
from .rateFactory import equilibrium, forwardRate, reverseRate


class ReactionSet(EntitySet):
    def __init__(self, mixture, reactions):
        EntitySet.__init__(self)

        from .Reaction import Reaction

        i = 0
        for reaction in reactions:
            i += 1
            K_c = equilibrium(mixture, reaction)
            k_f = forwardRate(mixture, reaction)
            k_r = reverseRate(mixture, reaction, K_c, k_f)
            r = Reaction(reaction, K_c, k_f, k_r)

            self.insert(i, r)

        return


# version
__id__ = "$Id$"

# End of file
