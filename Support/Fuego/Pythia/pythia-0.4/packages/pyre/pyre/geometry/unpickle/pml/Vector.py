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


class Vector(AbstractNode):

    tag = "vector"


    def content(self, content):
        self._vector = self._parse(content)
        return


    def notify(self, parent):
        parent.onVector(self._vector)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._vector = None
        return


# version
__id__ = "$Id$"

# End of file
