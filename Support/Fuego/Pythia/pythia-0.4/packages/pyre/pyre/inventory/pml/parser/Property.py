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


class Property(AbstractNode):


    tag = "property"


    def notify(self, parent):
        return parent.onProperty(self)


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self.name = attributes["name"]
        self.value = attributes.get("value")
        return
    

# version
__id__ = "$Id$"

# End of file 
