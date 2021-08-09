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


class Rotation(Transformation):


    def inspect(self, visitor):
        return visitor.onRotation(self)


    def __init__(self, body, vector, angle):
        Transformation.__init__(self, body)

        self.angle = angle
        self.vector = tuple(vector)

        self._info.log(str(self))

        return


    def __str__(self):
        return "rotation: body={%s}, vector=%r, angle=%s" % (self.body, self.vector, self.angle)


# version
__id__ = "$Id$"

#
# End of file
