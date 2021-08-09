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
from .Binary import Binary
import pyre.geometry.operations


class Difference(Binary):

    tag = "difference"


    def notify(self, parent):
        if not self._b1 or not self._b2:
            raise Exception("Not enough arguments for difference")
        
        difference = pyre.geometry.operations.subtract(self._b1, self._b2)
        parent.onDifference(difference)

        return


# version
__id__ = "$Id$"

# End of file
