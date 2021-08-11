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

import pyre.geometry.operations

from .Transformation import Transformation


class Translation(Transformation):

    tag = "translation"

    def onVector(self, vector):
        self._vector = vector
        return

    def notify(self, parent):
        if not self._body:
            raise Exception("Not enough arguments for translation")
        translation = pyre.geometry.operations.translate(
            body=self._body, vector=self._vector
        )
        parent.onTranslation(translation)
        return

    def __init__(self, root, attributes):
        Transformation.__init__(self, root, attributes)
        self._vector = None
        return


# version
__id__ = "$Id$"

# End of file
