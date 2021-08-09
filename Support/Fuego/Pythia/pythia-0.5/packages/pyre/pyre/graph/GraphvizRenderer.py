#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from builtins import object
class GraphvizRenderer(object):


    def draw(self, graph):
        self._graph = graph

        self._initialize()
        self._renderGraph()
        self._renderNodes()
        self._renderEdges()
        self._finalize()

        self._graph = None
        
        return self._text


    def __init__(self):
        self._graph = None
        self._text = ""
        return


    def _initialize(self):
        self._text = 'digraph "%s" {' % self._graph._id
        return


    def _finalize(self):
        self._text = self._text + "}"
        return


    def _renderGraph(self):
        if not self._graph._attributes: return ""
        self._text = self._text + " " + decorate(self._graph._attributes, ";") + ";"
        return


    def _renderNodes(self):
        for node in self._graph._nodes:

            if not node._id: text = " node"
            else: text = ' "%s"' % node._id
            
            if node._attributes: 
                text = text + " [" + decorate(node._attributes) + "]"
            text = text + ";"

            self._text = self._text + text

        return


    def _renderEdges(self):
        for edge in self._graph._edges:
            text = ' "%s" -> {"%s"}' % (edge._source.id(), edge._sink.id())
            if edge._attributes: 
                text = text + " [" + decorate(edge._attributes) + "]"
            text = text + ";"

            self._text = self._text + text

        return

# helpers

def decorate(table, delimiter=","):
    attributes = []

    for attribute, value in list(table.items()):
        attributes.append('%s="%s"' % (attribute, value))
    
    return join(attributes, delimiter)


# version
__id__ = "$Id$"

#
# End of file
