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

from AbstractNode import AbstractNode


class Composition(AbstractNode):


    def onBlock(self, body):
        self._setOperand(body)
        return


    def onCone(self, body):
        self._setOperand(body)
        return


    def onCylinder(self, body):
        self._setOperand(body)
        return


    def onPrism(self, body):
        self._setOperand(body)
        return


    def onPyramid(self, body):
        self._setOperand(body)
        return


    def onSphere(self, body):
        self._setOperand(body)
        return


    def onTorus(self, body):
        self._setOperand(body)
        return


    def onGeneralizedCone(self, body):
        self._setOperand(body)
        return


    def onDifference(self, body):
        self._setOperand(body)
        return


    def onIntersection(self, body):
        self._setOperand(body)
        return


    def onUnion(self, body):
        self._setOperand(body)
        return


    def onDilation(self, body):
        self._setOperand(body)
        return


    def onReflection(self, body):
        self._setOperand(body)
        return


    def onReversal(self, body):
        self._setOperand(body)
        return


    def onRotation(self, body):
        self._setOperand(body)
        return


    def onTranslation(self, body):
        self._setOperand(body)
        return


# version
__id__ = "$Id$"

# End of file
