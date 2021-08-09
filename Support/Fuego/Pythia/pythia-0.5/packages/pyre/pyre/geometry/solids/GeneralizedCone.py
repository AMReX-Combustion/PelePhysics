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


class GeneralizedCone(Primitive):
    def inspect(self, visitor):
        return visitor.onGeneralizedCone(self)

    def __init__(self, major, minor, scale, height):
        self.major = major
        self.minor = minor
        self.scale = scale
        self.height = height

        self._info.log("new %s" % self)

        return

    def __str__(self):
        return "cone: major=%s, minor=%s, scale=%s, height=%s" % (
            self.major,
            self.minor,
            self.scale,
            self.height,
        )


# version
__id__ = "$Id$"

#
# End of file
