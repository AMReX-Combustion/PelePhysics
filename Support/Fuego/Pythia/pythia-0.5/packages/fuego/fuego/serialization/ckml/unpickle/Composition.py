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
 

class Composition(AbstractNode):


    tag = "composition"


    def notify(self, parent):
        parent.onComposition(self._composition)
        return


    def onAtom(self, symbol, coefficient):
        self._composition.append((symbol, coefficient))
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._composition = []
        return
            

# version
__id__ = "$Id$"

#  End of file 
