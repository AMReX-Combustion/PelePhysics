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


from builtins import object


class Transport(object):
    def __init__(self, EPS, SIG, DIP, POL, ZROT, locator=None):
        self.eps = EPS
        self.sig = SIG
        self.dip = DIP
        self.pol = POL
        self.zrot = ZROT
        self.locator = locator
        return


# version
__id__ = "$Id$"

#  End of file
