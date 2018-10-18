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

from Binary import Binary
import pyre.geometry.operations


class Intersection(Binary):

    tag = "intersection"


    def notify(self, parent):
        if not self._b1 or not self._b2:
            raise "Not enough arguments for intersection"
        intersection = pyre.geometry.operations.intersect(self._b1, self._b2)
        parent.onIntersection(intersection)
        return


# version
__id__ = "$Id$"

# End of file
