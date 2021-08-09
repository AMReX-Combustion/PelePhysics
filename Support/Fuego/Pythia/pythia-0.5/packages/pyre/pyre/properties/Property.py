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


from pyre.inventory.Trait import Trait


class Property(Trait):


    def __init__(self, name, default=None, public=None, validator=None, tip="", doc=""):
        Trait.__init__(self, name, default=default, public=public, tip=tip, doc=doc)
        self.type = "property"
        self.validator = validator
        return


    def __set__(self, instance, value):
        value = self._cast(value)
        v = self.validator
        if v and not v(value):
            raise ValueError("value '%s' rejected by the '%s' validator" % (value, self.name))
            
        instance.__dict__[self.mangled] = value
        return


    def _cast(self, input):
        return input


# version
__id__ = "$Id$"

# End of file 
