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
 

class SRI(AbstractNode):


    tag = "sri"


    def notify(self, parent):
        parent.onSRI(self._parameters)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        a = float(attributes["a"])
        b = float(attributes["b"])
        c = float(attributes["c"])
        d = attributes.get("d")
        e = attributes.get("e")

        if d and e:
            self._parameters = (a, b, c, float(d), float(e))
        else:
            self._parameters = (a, b, c)
            
        return
            

# version
__id__ = "$Id$"

#  End of file 
