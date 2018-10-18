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

from Primitive import Primitive


class Torus(Primitive):


    def inspect(self, visitor):
        return visitor.onTorus(self)


    def __init__(self, major, minor):
        self.major = major
        self.minor = minor

        self._info.log("new %s" % self)
                 
        return
        

    def __str__(self):
        return "torus: major=%s, minor=%s" % (self.major, self.minor)

# version
__id__ = "$Id$"

#
# End of file
