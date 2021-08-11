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


class Sphere(AbstractNode):

    tag = "sphere"

    def notify(self, parent):
        sphere = pyre.geometry.solids.sphere(radius=self._radius)
        parent.onSphere(sphere)

        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._radius = self._parse(attributes["radius"])
        return


# version
__id__ = "$Id$"

# End of file
