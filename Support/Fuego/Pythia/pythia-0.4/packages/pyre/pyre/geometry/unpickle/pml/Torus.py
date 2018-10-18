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


class Torus(AbstractNode):

    tag = "torus"


    def notify(self, parent):
        torus = pyre.geometry.solids.torus(major=self._major, minor=self._minor)
        parent.onTorus(torus)

        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._major = self._parse(attributes["major"])
        self._minor = self._parse(attributes["minor"])
        return


# version
__id__ = "$Id$"

# End of file
