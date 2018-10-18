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

from pyre.components.Component import Component


class Device(Component):


    def createDevice(self):
        raise NotImplementedError("class '%s' must override 'device'" % self.__class__.__name__)


    def __init__(self, name):
        Component.__init__(self, name, "device")
        self.device = None
        return


    def _init(self, parent):
        device = self.createDevice()
        renderer = self.inventory.renderer.renderer
        device.renderer = renderer

        self.device = device
        
        return device


    class Inventory(Component.Inventory):

        from RendererFacility import RendererFacility

        inventory = (
            RendererFacility(),
            )


    


# version
__id__ = "$Id$"

# End of file 
