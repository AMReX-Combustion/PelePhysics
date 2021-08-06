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
import pyre.geometry.solids
from .AbstractNode import AbstractNode


class Pyramid(AbstractNode):

    tag = "pyramid"


    def notify(self, parent):
        pyramid = pyre.geometry.solids.pyramid()
        parent.onPyramid(pyramid)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        return


# version
__id__ = "$Id$"

# End of file
