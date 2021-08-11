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


def edgeAttributes():
    """return a list of valid attributes for Edge"""
    return list(Edge._validAttributes.keys())


class Edge(Drawable):
    def __init__(self, source, sink):
        Drawable.__init__(self)

        self._sink = sink
        self._source = source

        return

    _validAttributes = {
        "color": None,
        "decorate": None,
        "dir": None,
        "fontcolor": None,
        "fontname": None,
        "fontsize": None,
        "id": None,
        "label": None,
        "layer": None,
        "minlen": None,
        "style": None,
        "weight": None,
    }


# version
__id__ = "$Id$"

#
# End of file
