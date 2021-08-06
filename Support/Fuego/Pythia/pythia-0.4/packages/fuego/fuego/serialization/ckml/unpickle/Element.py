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
 

class Element(AbstractNode):


    tag = "element"


    def notify(self, parent):
        parent.onElement(self._element)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        symbol = attributes["id"]
        atomicWeight = attributes.get("atomicWeight")
        locator = self.documentNode().locator()

        self._element = self.documentNode().mechanism().newElement(symbol, atomicWeight, locator)
        return
            

# version
__id__ = "$Id$"

#  End of file 
