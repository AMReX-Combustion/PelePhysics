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

from Token import Token


class ReactionUnitsE(Token):

    
    pattern = "|".join(
        ["([Cc][Aa][Ll]/[Mm][Oo][Ll][Ee])",                   # cal/mole
         "([Kk][Cc][Aa][Ll]/[Mm][Oo][Ll][Ee])",               # kcal/mole
         "([Kk]?[Jj][Oo][Uu][Ll][Ee][Ss]/[Mm][Oo][Ll][Ee])",  # joule/mole
         "[Kk][Ee][Ll][Vv][Ii][Nn][Ss]",                      # kelvins
         "[Ee][Vv][Oo][Ll][Tt][Ss]"                           # ev
         ])


    def identify(self, auth):
        return auth.aReactionUnitsE(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.units_E = self.lexeme.lower()
        return


    def __str__(self):
         return "{Reaction units for E: %s}" % self.units_E


# version
__id__ = "$Id$"

#
# End of file
