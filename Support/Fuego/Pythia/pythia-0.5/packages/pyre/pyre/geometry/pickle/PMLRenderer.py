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

from pyre.geometry.Visitor import Visitor
from weaver.mills.XMLMill import XMLMill


class PMLRenderer(XMLMill, Visitor):


    # solids bodies
    def onBlock(self, block):
        line = '<block diagonal="(%s, %s, %s)"/>' % block.diagonal

        self._write(line)
        return


    def onCone(self, cone):
        line = '<cone height="%s" topRadius="%s" bottomRadius="%s"/>' % (
            cone.height, cone.top, cone.bottom)

        self._write(line)
        return


    def onCylinder(self, cylinder):
        line = '<cylinder height="%s" radius="%s"/>' % (
            cylinder.height, cylinder.radius)

        self._write(line)
        return


    def onPrism(self, prism):
        # NYI
        return self._abstract("onPrism")


    def onPyramid(self, pyramid):
        # NYI
        return self._abstract("onPyramid")


    def onSphere(self, sphere):
        line = '<sphere radius="%s"/>' % sphere.radius

        self._write(line)
        return


    def onTorus(self, torus):
        line = '<torus major="%s" minor="%s"/>' % (torus.major, torus.minor)

        self._write(line)
        return


    def onGeneralizedCone(self, cone):
        line = '<generalized-cone major="%s" minor="%s" scale="%s" height="%s"/>' % (
            cone.major, cone.minor, cone.scale, cone.height)

        self._write(line)
        return


    # Euler operations
    def onDifference(self, difference):
        self._write("<difference>")

        self._indent()
        difference.op1.inspect(self)
        difference.op2.inspect(self)
        self._outdent()

        self._write("</difference>")

        return


    def onIntersection(self, intersection):
        self._write("<intersection>")

        self._indent()
        intersection.op1.inspect(self)
        intersection.op2.inspect(self)
        self._outdent()

        self._write("</intersection>")

        return


    def onUnion(self, union):
        self._write("<union>")

        self._indent()
        union.op1.inspect(self)
        union.op2.inspect(self)
        self._outdent()

        self._write("</union>")

        return


    # transformations
    def onDilation(self, dilation):
        self._write("<dilation>")

        self._indent()
        body = dilation.body.inspect(self)
        self._write( "<scale>%g</scale>" % dilation.scale)
        self._outdent()

        self._write("</dilation>")
        return


    def onReflection(self, reflection):
        self._write("<reflection>")

        self._indent()
        body = reflection.body.inspect(self)
        self._write("<vector>(%s, %s, %s)</vector>" % reflection.vector)
        self._outdent()

        self._write("</reflection>")
        return


    def onReversal(self, reversal):
        self._write("<reversal>")

        self._indent()
        body = reversal.body.inspect(self)
        self._outdent()

        self._write("</reversal>")
        return


    def onRotation(self, rotation):
        self._write( "<rotation>")

        self._indent()
        rotation.body.inspect(self)
        self._write("<angle>%g</angle>" % rotation.angle)
        self._write("<vector>(%s, %s, %s)</vector>" % rotation.vector)
        self._outdent()

        self._write("</rotation>")
        return


    def onTranslation(self, translation):
        self._write("<translation>")

        self._indent()
        translation.body.inspect(self)
        self._write("<vector>(%s, %s, %s)</vector>" % translation.vector)
        self._outdent()

        self._write("</translation>")
        return


    def onGeometry(self, body):
        self._indent()
        body.inspect(self)
        self._outdent()
        return


    def __init__(self):
        XMLMill.__init__(self)
        return
            

    def _renderDocument(self, body, options=None):

        self._rep += ['', '<!DOCTYPE geometry>', '', '<geometry>', '' ]
        self.onGeometry(body)
        self._rep += ['</geometry>']

        return


# version
__id__ = "$Id$"

#
# End of file
