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

from pyre.inventory.Registry import Registry

from .AbstractNode import AbstractNode


class Weaver(AbstractNode):

    tag = "weaver"

    def onMill(self, mill):
        self._options.attach(mill)
        return

    def notify(self, parent):
        parent.onWeaver(self._options)
        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._options = Registry("root")

        return


# version
__id__ = "$Id$"

#  End of file