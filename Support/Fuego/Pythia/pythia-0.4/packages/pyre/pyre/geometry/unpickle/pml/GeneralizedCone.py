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


class GeneralizedCone(AbstractNode):

    tag = "generalized-cone"


    def notify(self, parent):
        cone = pyre.geometry.solids.generalizedCone(
            major=self._major, minor=self._minor, scale=self._scale,
            height=self._height)

        parent.onCone(cone)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._major = self._parse(attributes["major"])
        self._minor = self._parse(attributes["minor"])
        self._scale = self._parse(attributes["scale"])
        self._height = self._parse(attributes["height"])

        return


# version
__id__ = "$Id$"

# End of file
