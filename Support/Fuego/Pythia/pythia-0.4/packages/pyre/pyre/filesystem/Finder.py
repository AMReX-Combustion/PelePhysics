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
from .Inspector import Inspector


class Finder(Inspector):


    def find(self, root, regexp):
        import re
        self._target = re.compile(regexp)
        self._nodes = []

        root.id(self)

        return self._nodes


    def onDirectory(self, node):
        for entry in node.children():
            if self._target.match(entry.name):
                self._nodes.append(entry)

        for entry in node.subdirectories():
            entry.id(self)

        return


    def __init__(self):
        self._nodes = []
        self._target = None
        return


# version
__id__ = "$Id$"

#  End of file 
