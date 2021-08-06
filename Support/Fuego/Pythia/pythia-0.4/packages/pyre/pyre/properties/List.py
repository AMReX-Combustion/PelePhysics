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


class List(Property):


    def __init__(self, name, default=[], public=None, validator=None, tip="", doc=""):
        super(List, self).__init__(name, default, public, validator, tip, doc)
        self.type = "list"
        return


    def _cast(self, value):
        if isinstance(value, str):
            try:
                value = eval(value)
            except:
                raise TypeError("property '%s': could not convert '%s' to a list" % (
                    self.name, value))

        if isinstance(value, list):
            return value
            
        raise TypeError("property '%s': could not convert '%s' to a list" % (self.name, value))
    

# version
__id__ = "$Id$"

# End of file 
