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


class Cylinder(Primitive):
    def inspect(self, visitor):
        return visitor.onCylinder(self)

    def __init__(self, radius, height):
        self.radius = radius
        self.height = height

        self._info.log("new %s" % self)

        return

    def __str__(self):
        return "cylinder: radius=%s, height=%s" % (self.radius, self.height)


# version
__id__ = "$Id$"

#
# End of file
