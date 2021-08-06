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

from __future__ import absolute_import
from .Drawable import Drawable
from .Node import Node
from .Edge import Edge


class Graph(Drawable):


    # nameless nodes set defaults for subsequent nodes
    def node(self, id=""):
        node = Node(id)
        self._nodes.append(node)
        return node


    def connect(self, source, sinks):
        edge = Edge(source, sinks)
        self._edges.append(edge)
        return edge


    def __init__(self, name):
        Drawable.__init__(self)

        self._graph = ""
        self._nodes = []
        self._edges = []
        self._id = name

        return


    _validAttributes = {
        "center" : ["true", "false"],
        "clusterrank" : ["local", "global", "none"],
        "color" : None,
        "concentrate": ["true", "false"],
        "fontcolor" : None,
        "fontname" : None,
        "fontsize" : None,
        "label" : None,
        "layers" : None,
        "margin" : None,
        "mclimit" : None,
        "nodesep" : None,
        "nslimit" : None,
        "ordering" : None,
        "orientation" : None,
        "page" : None,
        "rank" : None,
        "rankdir" : None,
        "ranksep" : None,
        "ratio" : None,
        "size" : None
        }


# version
__id__ = "$Id$"

#
# End of file
