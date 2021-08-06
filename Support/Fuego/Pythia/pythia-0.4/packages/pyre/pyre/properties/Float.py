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
from .Property import Property


class Float(Property):


    def __init__(self, name, default=0.0, public=None, validator=None, tip="", doc=""):
        super(Float, self).__init__(name, default, public, validator, tip, doc)
        self.type = "float"
        return


    def _cast(self, value):
        return float(value)
    

# version
__id__ = "$Id$"

# End of file 
