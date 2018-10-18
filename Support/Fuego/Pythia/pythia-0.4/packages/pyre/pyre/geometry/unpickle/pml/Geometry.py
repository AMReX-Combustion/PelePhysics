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

from Composition import Composition


class Geometry(Composition):

    tag = "geometry"


    def notify(self, parent):
        parent.onGeometry(self._bodies)
        return


    def __init__(self, root, attributes):
        Composition.__init__(self, root, attributes)
        self._bodies = []
        return


    def _setOperand(self, body):
        self._bodies.append(body)
        return


# version
__id__ = "$Id$"

# End of file
