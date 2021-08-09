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

from pyre.units import length
from pyre.xml.Node import Node


class AbstractNode(Node):
    def __init__(self, root, attributes):
        Node.__init__(self, root)
        return

    def _parse(self, expr):
        return eval(expr, self._lengthUnits)

    _lengthUnits = length.__dict__


# version
__id__ = "$Id$"

#
# End of file
