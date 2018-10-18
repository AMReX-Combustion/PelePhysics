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

import re

from Token import Token
from RegularExpressions import eol, whitespaceOpt, species, inlineNumber


class ReactionEfficiencies(Token):

    # don't be tempted to optimize prematurely and swallow the eol as well:
    #     comments may follow the paramaters on this line
    pattern = r"(" \
             + species + whitespaceOpt + inlineNumber + whitespaceOpt \
             + r")+" #+ eol

    def identify(self, auth):
        return auth.someReactionEfficiencies(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)

        self.efficiencies = []
        text = self.lexeme

        i = 0
        split = text.split("/")
        while i < len(split) - 1:
            species = split[i].strip()
            try:
                efficiency = float(split[i+1])
            except ValueError:
                # this can't happen because the regexp requires a float here 
                import pyre
                str = "Could not convert '%s' into a number" % split[i+1]
                pyre.debug.Firewall.hit(str)
                return

            self.efficiencies.append( (species, efficiency) )
            i = i + 2
        
        return


    def __str__(self):
        return "{third body efficiencies}"


# version
__id__ = "$Id$"

#
# End of file
