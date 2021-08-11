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


class ReactionTROE(Token):

    pattern = r"[Tt][Rr][Oo][Ee]"

    def identify(self, auth):
        return auth.aReactionTROE(self)

    def __str__(self):
        return "{TROE fall off}"


# version
__id__ = "$Id$"

#
# End of file
