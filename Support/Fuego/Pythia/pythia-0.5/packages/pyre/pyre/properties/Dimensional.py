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
from builtins import range
from .Property import Property


class Dimensional(Property):


    def __init__(self, name, default=0.0, public=None, validator=None, tip="", doc=""):
        super(Dimensional, self).__init__(name, default, public, validator, tip, doc)
        self.type = "dimensional"

        try:
            self.len = len(default)
        except TypeError:
            self.len = 0
            
        return


    def _cast(self, value):
        candidate = value
        if isinstance(value, str):
            import pyre.units
            parser = pyre.units.parser()
            candidate = parser.parse(value)

        self._checkDimensions(candidate, value)

        return value


    def _checkDimensions(self, candidate, setting):
        try:
            size = len(candidate)
        except TypeError:
            size = 0
        
        if size != self.len:
            raise ValueError("value '%s' is not the same shape as the default '%s'" % (
                setting, self.default))
        
        if self.len == 0:
            tokens = [candidate]
            target = [self.default]
        else:
            tokens = candidate
            target = self.default

        for i in range(size):
            if tokens[i].derivation != target[i].derivation:
                raise ValueError("dimension mismatch between input '%s' and target '%s'" % (
                    setting, self.default))

        return


# version
__id__ = "$Id$"

# End of file 
