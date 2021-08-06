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
 

class LowPressure(AbstractNode):


    tag = "low"


    def notify(self, parent):
        parent.onLowPressure(self._parameters)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        A = float(attributes["A"])
        beta = float(attributes["beta"])
        E = float(attributes["E"])

        self._parameters = (A, beta, E)
        return
            

# version
__id__ = "$Id$"

#  End of file 
