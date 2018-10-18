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


import journal
from pyre.parsing.Locator import Locator
from Node import Node


class DocumentNode(Node):


    def source(self):
        return self._source


    def locator(self):
        file = self._source
        line = self._locator.getLineNumber()
        column = self._locator.getColumnNumber()

        return Locator(file, line, column)


    # the node that represents the root of the document
    def documentNode(self):
        return self


    def node(self, name, attributes):
        nodeFactory = self._dtd.get(name)
        if not nodeFactory:
            firewall = journal.firewall("pyre.xml.parsing")
            firewall.log("tag '%s' has no registered handler" % name)
            return

        return nodeFactory(self, attributes)


    def setLocator(self, locator):
        self._locator = locator
        return


    def __init__(self, source, nodes):
        Node.__init__(self, None)

        self._source = source
        self._locator = None
        self._dtd = createDispatchTable(nodes)

        return


# helpers

def createDispatchTable(nodes):

    dispatcher = {}
    for node in nodes:
        dispatcher[node.tag] = node

    return dispatcher



# version
__id__ = "$Id$"

#  End of file 
