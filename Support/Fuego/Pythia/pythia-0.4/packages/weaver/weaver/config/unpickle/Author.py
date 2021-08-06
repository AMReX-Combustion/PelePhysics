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
from .AbstractNode import AbstractNode


class Author(AbstractNode):


    tag = "author"


    def content(self, text):
        self._author = text
        return


    def notify(self, parent):
        parent.onAuthor(self._author)
        return


    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._author = ""
        return


# version
__id__ = "$Id$"

#  End of file 
