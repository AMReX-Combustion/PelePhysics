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
from .Token import Token


class ReactionFORD(Token):


    pattern = r"[Ff][Oo][Rr][Dd]"


    def identify(self, auth):
        return auth.aReactionFORD(self)


    def __str__(self): 
        return "{Reaction order}"


# version
__id__ = "$Id$"

#
# End of file
