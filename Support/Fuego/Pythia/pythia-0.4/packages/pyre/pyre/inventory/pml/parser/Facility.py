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


class Facility(AbstractNode):


    tag = "facility"


    def notify(self, parent):
        return parent.onFacility(self.facility)


    def onFacility(self, facility):
        self.facility.attach(facility)
        return


    def onProperty(self, property):
        self.facility.set(property.name, property.value)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        name = attributes["name"]

        from pyre.inventory.Registry import Registry
        self.facility = Registry(name)
        return
    

# version
__id__ = "$Id$"

# End of file 
