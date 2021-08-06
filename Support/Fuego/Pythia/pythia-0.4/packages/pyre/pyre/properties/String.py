#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import
from builtins import str
from .Property import Property


class String(Property):


    def __init__(self, name, default="", public=None, validator=None, tip="", doc=""):
        super(String, self).__init__(name, default, public, validator, tip, doc)
        self.type = "str"
        return


    def _cast(self, value):
        return str(value)


# version
__id__ = "$Id$"

# End of file 
