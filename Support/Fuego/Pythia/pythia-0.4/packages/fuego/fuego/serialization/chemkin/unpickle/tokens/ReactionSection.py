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


class ReactionSection(Token):


    pattern = r"[Rr][Ee][Aa][Cc]([Tt][Ii][Oo][Nn][Ss])?"


    def identify(self, auth):
        return auth.aReactionSection(self)


    def __str__(self):
         return "{Reaction section}"


# version
__id__ = "$Id$"

#
# End of file
