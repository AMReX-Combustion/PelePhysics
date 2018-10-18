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
from RegularExpressions import species


class SpeciesName(Token):


    pattern = species


    def identify(self, auth):
        return auth.aSpeciesName(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.name = self.lexeme
        return

    
    def __str__(self):
         return "{species name=" + self.name + "}"


# version
__id__ = "$Id$"

#
# End of file
