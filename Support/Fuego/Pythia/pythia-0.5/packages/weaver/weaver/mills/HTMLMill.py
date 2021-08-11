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

from .BlockMill import BlockMill


class HTMLMill(BlockMill):

    names = ["html"]

    def __init__(self):
        BlockMill.__init__(
            self,
            "<!--",
            " !",
            "-->",
            '<!doctype html public "-//w3c//dtd html 4.0 transitional//en">',
        )
        return


# version
__id__ = "$Id$"

#  End of file
