#!/usr/bin/env python
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
#  <LicenseText>
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import

from .Primitive import Primitive


class Block(Primitive):
    def inspect(self, visitor):
        return visitor.onBlock(self)

    def __init__(self, diagonal):
        self.diagonal = tuple(diagonal)

        self._info.log("new %s" % self)

        return

    def __str__(self):
        return "block: diagonal=(%s, %s, %s)" % self.diagonal


# version
__id__ = "$Id$"

#
# End of file
