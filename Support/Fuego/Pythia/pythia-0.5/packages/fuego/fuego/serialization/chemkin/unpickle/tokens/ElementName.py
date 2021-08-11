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

from .RegularExpressions import namedElement, namedInlineNumber, whitespaceOpt
from .Token import Token


class ElementName(Token):

    _patternNames = ("element_name", "element_weight")
    _patternTemplate = (
        namedElement + whitespaceOpt + r"(" + namedInlineNumber + r")?"
    )
    pattern = _patternTemplate % _patternNames

    def identify(self, auth):
        return auth.anElementName(self)

    def __init__(self, match, groups):
        Token.__init__(self, match, groups)
        self.name = groups["element_name"]
        self.weight = groups["element_weight"]
        return

    def __str__(self):
        str = "{element: name=" + self.name
        if self.weight:
            str = str + ", weight=" + self.weight
        str = str + "}"
        return str


# version
__id__ = "$Id$"

#
# End of file
