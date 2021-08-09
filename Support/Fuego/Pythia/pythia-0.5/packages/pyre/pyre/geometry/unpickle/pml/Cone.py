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


class Cone(AbstractNode):

    tag = "cone"

    def notify(self, parent):
        cone = pyre.geometry.solids.cone(
            top=self._topRadius, bottom=self._bottomRadius, height=self._height
        )

        parent.onCone(cone)

        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._topRadius = self._parse(attributes["topRadius"])
        self._bottomRadius = self._parse(attributes["bottomRadius"])

        self._height = self._parse(attributes["height"])

        return


# version
__id__ = "$Id$"

# End of file
