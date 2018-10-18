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

from Drawable import Drawable


def nodeAttributes():
    """return a list of valid attributes for Node"""
    return Node._validAttributes.keys()


class Node(Drawable):


    def id(self): return self._id


    def __init__(self, id):
        Drawable.__init__(self)
        self._id = id
        return


    _validAttributes = {
        "color" : None,
        "fontcolor" : None,
        "fontname" : None,
        "fontsize" : None,
        "height" : None,
        "label" : None,
        "layer" : None,
        "shape" : None,
        "shapefile" : None,
        "style" : None,
        "width" : None
        }

# version
__id__ = "$Id$"

#
# End of file
