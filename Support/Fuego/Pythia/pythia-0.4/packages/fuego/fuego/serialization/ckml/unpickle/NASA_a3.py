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
 

class NASA_a3(AbstractNode):


    tag = "a3"


    def notify(self, parent):
        parent.onNASA_a3(self._value)
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
