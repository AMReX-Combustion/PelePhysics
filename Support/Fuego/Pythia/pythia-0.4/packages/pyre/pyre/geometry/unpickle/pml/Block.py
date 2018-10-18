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


class Block(AbstractNode):

    tag = "block"


    def notify(self, parent):
        block = pyre.geometry.solids.block(diagonal=self._diagonal)
        parent.onBlock(block)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._diagonal = self._parse(attributes["diagonal"])
        return


# version
__id__ = "$Id$"

# End of file
