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


class ReactionDuplicate(Token):
    
    pattern = r"[Dd][Uu][Pp]([Ll][Ii][Cc][Aa][Tt][Ee])?"


    def identify(self, auth):
        return auth.aReactionDuplicate(self)


    def __str__(self):
        return "{duplicate reaction}"


# version
__id__ = "$Id$"

#
# End of file
