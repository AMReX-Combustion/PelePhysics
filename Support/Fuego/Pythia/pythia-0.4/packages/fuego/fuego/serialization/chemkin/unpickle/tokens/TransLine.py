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
from RegularExpressions import eol, whitespaceOpt


class TransLine(Token):


    pattern = r"(?P<tran_line>.{19,19})(?P<tran_type>[0-2])(?P<tran_line_2>.*$)"


    def identify(self, auth):
        return auth.aTransLine(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.text = groups["tran_line"]
        self.id = int(groups["tran_type"])
        self.text_2 = groups["tran_line_2"]

        return


    def __str__(self):
        return "{trans line: %s}" % self.text.strip()


# version
__id__ = "$Id$"

#
# End of file
