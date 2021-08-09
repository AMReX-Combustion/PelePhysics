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

from .LineComments import LineComments
from .Mill import Mill


class LineMill(Mill, LineComments):
    def __init__(self, comment, firstLine):
        Mill.__init__(self)
        LineComments.__init__(self)

        self.commentLine = comment
        self.firstLine = firstLine

        return


# version
__id__ = "$Id$"

#  End of file
