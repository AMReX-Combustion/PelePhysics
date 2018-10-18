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
 

class ReactionUnits(AbstractNode):


    tag = "units"


    def notify(self, parent):
        parent.onReactionUnits(self._reactionUnits)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)

        self._reactionUnits = {
            "activation": attributes.get("activation", "cal/mole"),
            "prefactor": attributes.get("prefactor", "mole/cm**3")
            }

        return
            

# version
__id__ = "$Id$"

#  End of file 
