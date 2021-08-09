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

from .Transport import Transport


class TRANLIN(Transport):
    def __init__(self, EPS, SIG, DIP, POL, ZROT, locator=None):
        Transport.__init__(self, EPS, SIG, DIP, POL, ZROT, locator)
        self.parameters = []
        return


# version
__id__ = "$Id$"

#  End of file
