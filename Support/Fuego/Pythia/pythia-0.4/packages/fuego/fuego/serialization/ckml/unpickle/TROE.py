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
 

class TROE(AbstractNode):


    tag = "troe"


    def notify(self, parent):
        parent.onTROE(self._parameters)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        a = float(attributes["a"])
        Ts = float(attributes["Ts"])
        T3s = float(attributes["T3s"])

        T2s = attributes.get("T2s")

        if T2s:
            self._parameters = (a, T3s, Ts, float(T2s))
        else:
            self._parameters = (a, T3s, Ts)
            
        return
            

# version
__id__ = "$Id$"

#  End of file 
