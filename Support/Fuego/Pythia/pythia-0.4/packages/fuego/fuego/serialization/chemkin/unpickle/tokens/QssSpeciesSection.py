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


class QssSpeciesSection(Token):


    pattern = r"[Qq][Ss][Ss]?"


    def identify(self, auth):
        return auth.aQssSpeciesSection(self)


    def __str__(self): 
        return "{QSS Species section}"


# version
__id__ = "$Id$"

#
# End of file
