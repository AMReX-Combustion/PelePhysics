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


from AbstractNode import AbstractNode


class Component(AbstractNode):


    tag = "component"


    def notify(self, parent):
        return parent.onComponent(self.component)


    def onComponent(self, component):
        self.component.attach(component)
        return


    def onFacility(self, facility):
        self.facility.attach(facility)
        return


    def onProperty(self, property):
        self.component.set(property.name, property.value)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        name = attributes["name"]

        from pyre.inventory.Registry import Registry
        self.component = Registry(name)
        return
    

# version
__id__ = "$Id$"

# End of file 
