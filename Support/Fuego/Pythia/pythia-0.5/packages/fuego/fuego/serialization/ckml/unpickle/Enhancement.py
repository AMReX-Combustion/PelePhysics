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


class Enhancement(AbstractNode):

    tag = "enhancement"

    def notify(self, parent):
        parent.onEnhancement(self._species, self._coefficient)
        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._species = attributes["species"]
        self._coefficient = float(attributes["factor"])

        return


# version
__id__ = "$Id$"

#  End of file
