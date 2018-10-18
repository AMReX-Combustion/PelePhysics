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

from pyre.xml.DocumentNode import DocumentNode


class Document(DocumentNode):


    def bodies(self):
        return self._bodies


    def onGeometry(self, bodies):
        self._bodies = bodies
        return


    def __init__(self, source):
        DocumentNode.__init__(self, source, documentNodes())

        self._bodies = []

        return


# helpers

def documentNodes():

    from Geometry import Geometry

    # solids
    from Block import Block
    from Cone import Cone
    from Cylinder import Cylinder
    from Prism import Prism
    from Pyramid import Pyramid
    from Sphere import Sphere
    from Torus import Torus
    from GeneralizedCone import GeneralizedCone

    # compositions
    from Difference import Difference
    from Intersection import Intersection
    from Union import Union

    # transformations
    from Dilation import Dilation
    from Reflection import Reflection
    from Reversal import Reversal
    from Rotation import Rotation
    from Translation import Translation

    # data
    from Angle import Angle
    from Scale import Scale
    from Vector import Vector
    
    nodes = [
        Geometry,

        # solids
        Block, Cone, Cylinder, Prism, Pyramid, Sphere, Torus,
        GeneralizedCone,

        # compositions
        Difference, Intersection, Union,

        # transformations
        Dilation, Reflection, Reversal, Rotation, Translation,

        # data
        Angle, Scale, Vector,

        ]

    return nodes


# version
__id__ = "$Id$"

# End of file
