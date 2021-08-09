#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from __future__ import absolute_import
from .AbstractNode import AbstractNode


class Registry(AbstractNode):


    tag = "registry"


    def notify(self, parent):
        parent.onRegistry(self.registry)
        return


    def onComponent(self, component):
        self.registry.attach(component)
        return


    def onFacility(self, facility):
        self.registry.attach(facility)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        from pyre.inventory.Registry import Registry
        self.component = Registry("root")

        return


# version
__id__ = "$Id$"

# End of file 
