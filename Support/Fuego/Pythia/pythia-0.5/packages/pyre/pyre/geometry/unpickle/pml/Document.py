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

    # data
    from .Angle import Angle

    # solids
    from .Block import Block
    from .Cone import Cone
    from .Cylinder import Cylinder

    # compositions
    from .Difference import Difference

    # transformations
    from .Dilation import Dilation
    from .GeneralizedCone import GeneralizedCone
    from .Geometry import Geometry
    from .Intersection import Intersection
    from .Prism import Prism
    from .Pyramid import Pyramid
    from .Reflection import Reflection
    from .Reversal import Reversal
    from .Rotation import Rotation
    from .Scale import Scale
    from .Sphere import Sphere
    from .Torus import Torus
    from .Translation import Translation
    from .Union import Union
    from .Vector import Vector

    nodes = [
        Geometry,
        # solids
        Block,
        Cone,
        Cylinder,
        Prism,
        Pyramid,
        Sphere,
        Torus,
        GeneralizedCone,
        # compositions
        Difference,
        Intersection,
        Union,
        # transformations
        Dilation,
        Reflection,
        Reversal,
        Rotation,
        Translation,
        # data
        Angle,
        Scale,
        Vector,
    ]

    return nodes


# version
__id__ = "$Id$"

# End of file
