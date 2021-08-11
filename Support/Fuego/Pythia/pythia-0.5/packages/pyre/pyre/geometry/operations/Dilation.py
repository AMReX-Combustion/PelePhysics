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

from builtins import str

from .Transformation import Transformation


class Dilation(Transformation):
    def inspect(self, visitor):
        return visitor.onDilation(self)

    def __init__(self, body, scale):
        Transformation.__init__(self, body)
        self.scale = scale

        self._info.log(str(self))

        return

    def __str__(self):
        return "dilation: body={%s}, scale=%s" % (self.body, self.scale)


# version
__id__ = "$Id$"

#
# End of file
