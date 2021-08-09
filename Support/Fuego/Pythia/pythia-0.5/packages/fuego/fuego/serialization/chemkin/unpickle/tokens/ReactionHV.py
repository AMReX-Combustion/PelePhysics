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

from .RegularExpressions import namedInlineNumber, whitespaceOpt
from .Token import Token


class ReactionHV(Token):

    pattern = r"[Hh][Vv]"

    def identify(self, auth):
        return auth.aReactionHV(self)

    def __str__(self):
        return "{Radiation}"


# version
__id__ = "$Id$"

#
# End of file
