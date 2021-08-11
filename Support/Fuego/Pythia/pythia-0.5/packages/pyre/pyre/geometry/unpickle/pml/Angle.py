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

from .AbstractNode import AbstractNode


class Angle(AbstractNode):

    tag = "angle"

    def content(self, content):
        self._angle = float(content)
        return

    def notify(self, parent):
        parent.onAngle(self._angle)
        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._angle = None
        return


# version
__id__ = "$Id$"

# End of file
