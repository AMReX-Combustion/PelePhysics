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

import pyre.geometry.solids
from AbstractNode import AbstractNode


class Prism(AbstractNode):
    
    tag = "prism"


    def notify(self, parent):
        prism = pyre.geometry.solids.prism()
        parent.onPrism(prism)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        return


# version
__id__ = "$Id$"

# End of file
