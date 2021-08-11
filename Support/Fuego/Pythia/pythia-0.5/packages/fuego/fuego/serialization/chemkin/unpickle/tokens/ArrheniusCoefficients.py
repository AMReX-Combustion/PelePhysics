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

#
# It appears that this token class does not get used
# so I placed a Firewall in the constructor
#

from __future__ import absolute_import

from builtins import map

from .RegularExpressions import namedNumbers_3
from .Token import Token


class ArrheniusCoefficients(Token):

    pattern = namedNumbers_3 % ("arrhenius_A", "arrhenius_beta", "arrhenius_E")

    def identify(self, auth):
        return auth.someArrheniusCoefficients(self)

    def __init__(self, match, groups):
        Token.__init__(self, match, groups)

        import pyre

        pyre.debug.Firewall.hit("arrhenius coefficients token")

        numbers = list(map(groups.get, arrhenius))

        try:
            self.parameters = list(map(float, numbers))
        except ValueError:
            # this can't happen because the regexp requires three floats
            import pyre

            msg = "Could not convert /%s/ into a list of numbers" % text
            pyre.debug.Firewall.hit(msg)
            return

        return

    def __str__(self):
        return "{Arrhenius coefficients: %s}" % self.parameters


# version
__id__ = "$Id$"

#
# End of file
