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
import pyre.geometry.solids
from .AbstractNode import AbstractNode


class Cylinder(AbstractNode):

    tag = "cylinder"


    def notify(self, parent):
        cylinder = pyre.geometry.solids.cylinder(radius=self._radius, height=self._height)
        parent.onCylinder(cylinder)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._radius = self._parse(attributes["radius"])
        self._height = self._parse(attributes["height"])
        return


# version
__id__ = "$Id$"

# End of file
