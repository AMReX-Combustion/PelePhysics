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

import pyre.geometry.operations
from Transformation import Transformation


class Rotation(Transformation):

    tag = "rotation"


    def onAngle(self, angle):
        self._angle = angle
        return


    def onVector(self, vector):
        self._vector = vector
        return


    def notify(self, parent):
        if not self._body:
            raise "Not enough arguments for rotation"
        rotation = pyre.geometry.operations.rotate(
            body=self._body, vector=self._vector, angle=self._angle)
        parent.onRotation(rotation)
        return


    def __init__(self, root, attributes):
        Transformation.__init__(self, root, attributes)
        self._angle = None
        self._vector = None
        return


# version
__id__ = "$Id$"

# End of file
