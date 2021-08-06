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
from builtins import map
from .Token import Token
from .RegularExpressions import eol, whitespace, whitespaceOpt, namedNumbers_3


class TemperatureRange(Token):

    _patternNames = ("lowT", "commonT", "highT")

    pattern = whitespaceOpt \
             + namedNumbers_3 % _patternNames \
             + whitespaceOpt + eol


    def identify(self, auth):
        return auth.aTemperatureRange(self)


    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        text = list(map(groups.get, self._patternNames))
        try:
            self.range = tuple(map(float, text))
        except ValueError:
            # this can't happen because the regexp requires a float here 
            import pyre
            msg = "Can't convert '%s' into a list of numbers" % text
            pyre.debug.Firewall.hit(msg)
            return

        return


    def __str__(self):
        return "{temperature range: <%g><%g><%g>}" % (self.range[0], self.range[1], self.range[2])
    

# version
__id__ = "$Id$"

#
# End of file
