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

from .RegularExpressions import whitespace
from .Token import Token


class ThermoSection(Token):

    pattern = (
        r"[Tt][Hh][Ee][Rr]([Mm][Oo])?("
        + whitespace
        + "(?P<thermo_all>[Aa][Ll][Ll]))?"
    )

    def identify(self, auth):
        return auth.aThermoSection(self)

    def __init__(self, match, groups):
        Token.__init__(self, match, groups)

        if groups["thermo_all"]:
            self._all = 1
        else:
            self._all = 0
        return

    def __str__(self):
        return "{Thermo section: all=%d}" % self._all


# version
__id__ = "$Id$"

#
# End of file
