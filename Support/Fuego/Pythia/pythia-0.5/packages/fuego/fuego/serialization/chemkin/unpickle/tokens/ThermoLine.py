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

from .RegularExpressions import eol, whitespaceOpt
from .Token import Token


class ThermoLine(Token):

    pattern = (
        r"(?P<thermo_line>.{79,79})(?P<thermo_type>[1-4])"
        + whitespaceOpt
        + eol
    )

    def identify(self, auth):
        return auth.aThermoLine(self)

    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.text = groups["thermo_line"]
        self.id = int(groups["thermo_type"])

        return

    def __str__(self):
        return "{thermo line: %d}" % self.id


# version
__id__ = "$Id$"

#
# End of file
