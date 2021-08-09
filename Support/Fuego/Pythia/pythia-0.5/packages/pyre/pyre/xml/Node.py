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

from builtins import object
import journal


class Node(object):

    # descendants must specify their tag by setting self.tag
    tag = "unknown"


    # the node that represents the root of the document
    def documentNode(self):
        return self._documentNode


    # abstract methods
    def content(self, text):
        raise NotImplementedError(
            "class '%s' should override method 'content'" % self.__class__.__name__)


    def notify(self, target):
        raise NotImplementedError(
            "class '%s' should override method 'notify'" % self.__class__.__name__)


    def __init__(self, documentNode):
        self._documentNode = documentNode
        return


    _info = journal.debug("pyre.xml.parsing")


# version
__id__ = "$Id$"

#  End of file 
