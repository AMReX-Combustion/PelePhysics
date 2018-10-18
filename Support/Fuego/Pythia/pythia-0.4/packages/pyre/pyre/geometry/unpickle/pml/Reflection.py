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


class Reflection(Transformation):

    tag = "reflection"


    def onVector(self, vector):
        self._vector = vector
        return


    def notify(self, parent):
        if not self._body:
            raise "Not enough arguments for reflection"
        reflection = pyre.geometry.operations.reflect(body=self._body, vector=self._vector)
        parent.onReflection(reflection)
        return


    def __init__(self, root, attributes):
        Transformation.__init__(self, root, attributes)
        self._vector = None
        return


# version
__id__ = "$Id$"

# End of file
