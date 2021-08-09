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


class SpeciesSection(Token):

    pattern = r"[Ss][Pp][Ee][Cc]([Ii][Ee][Ss])?"

    def identify(self, auth):
        return auth.aSpeciesSection(self)

    def __str__(self):
        return "{Species section}"


# version
__id__ = "$Id$"

#
# End of file
