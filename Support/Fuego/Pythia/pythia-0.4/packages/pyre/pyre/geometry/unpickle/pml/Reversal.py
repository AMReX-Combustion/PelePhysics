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
import pyre.geometry.operations
from .Transformation import Transformation


class Reversal(Transformation):

    tag = "reversal"


    def notify(self, parent):
        if not self._body:
            raise Exception("Not enough arguments for reversal")
        reversal = pyre.geometry.operations.reverse(body=self._body)
        parent.onReversal(reversal)
        return


# version
__id__ = "$Id$"

# End of file
