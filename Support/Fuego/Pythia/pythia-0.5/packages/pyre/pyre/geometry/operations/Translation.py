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

from builtins import str

from .Transformation import Transformation


class Translation(Transformation):
    def inspect(self, visitor):
        return visitor.onTranslation(self)

    def __init__(self, body, vector):
        Transformation.__init__(self, body)
        self.vector = tuple(vector)

        self._info.log(str(self))

        return

    def __str__(self):
        return "translation: body={%s}, vector=%r" % (self.body, self.vector)


# version
__id__ = "$Id$"

#
# End of file
