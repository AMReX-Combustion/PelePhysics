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
 

class Reagents(AbstractNode):


    tag = "reagents"


    def notify(self, parent):
        parent.onReagents(self._reactants, self._products)
        return


    def onReactant(self, species, coefficient):
        self._reactants.append((species, coefficient))
        return


    def onProduct(self, species, coefficient):
        self._products.append((species, coefficient))
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._reactants = []
        self._products = []
        return
            

# version
__id__ = "$Id$"

#  End of file 
