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

from AbstractNode import AbstractNode
 

class NASA_a1(AbstractNode):


    tag = "a1"


    def notify(self, parent):
        parent.onNASA_a1(self._value)
        return


    def content(self, text):
        self._value = float(text)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._value = 0.0
        return
            

# version
__id__ = "$Id$"

#  End of file 
