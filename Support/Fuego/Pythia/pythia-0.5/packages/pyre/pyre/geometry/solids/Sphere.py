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


class Sphere(Primitive):


    def inspect(self, visitor):
        return visitor.onSphere(self)


    def __init__(self, radius):
        self.radius = radius

        self._info.log("new %s" % self)
                 
        return


    def __str__(self):
        return "sphere: radius=%s" % self.radius


# version
__id__ = "$Id$"

#
# End of file
