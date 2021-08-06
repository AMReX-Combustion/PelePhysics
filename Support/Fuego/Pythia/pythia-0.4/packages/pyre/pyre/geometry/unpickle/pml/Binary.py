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
from .Composition import Composition


class Binary(Composition):


    def __init__(self, root, attributes=None):
        Composition.__init__(self, root, attributes)
        
        self._b1 = None
        self._b2 = None

        return


    def _setOperand(self, body):
        if not self._b1:
            self._b1 = body
        elif not self._b2:
            self._b2 = body
        else:
            raise "Too many nested tags in binary body"
            
        return


# version
__id__ = "$Id$"

# End of file
