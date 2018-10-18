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
from RegularExpressions import whitespace


class TransSection(Token):


    pattern = r"[Tt][Rr][Aa][Nn]([Ss])?(" + whitespace + "(?P<trans_all>[Aa][Ll][Ll]))?"


    def identify(self, auth): 
        return auth.aTransSection(self)


    #def __init__(self, match, groups):
    #    Token.__init__(self, match, groups)

    #    if groups["trans_all"]:
    #        self._all = 1
    #    else:
    #        self._all = 0
    #    return


    def __str__(self):
        return "{Trans section}"

    
# version
__id__ = "$Id$"

#
# End of file
