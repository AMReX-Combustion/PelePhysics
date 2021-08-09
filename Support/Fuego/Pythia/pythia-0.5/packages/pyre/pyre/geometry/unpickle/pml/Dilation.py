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


class Dilation(Transformation):

    tag = "dilation"


    def onScale(self, scale):
        self._scale = scale
        return


    def notify(self, parent):
        if not self._body:
            raise Exception("Not enough arguments for dilation")
        dilation = pyre.geometry.operations.dilate(body=self._body, scale=self._scale)
        parent.onDilation(dilation)
        return


    def __init__(self, root, attributes):
        Transformation.__init__(self, root, attributes)
        self._scale = None
        return


# version
__id__ = "$Id$"

# End of file
