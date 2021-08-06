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
from builtins import object
class Device(object):


    def record(self, entry):
        str = self.renderer.render(entry)
        self._write(str)
        return


    def __init__(self, renderer=None):
        if renderer is None:
            from .Renderer import Renderer
            renderer = Renderer()

        self.renderer = renderer

        return


    def _write(self, entry):
        raise NotImplementedError("class '%s' must override '_write'" % self.__class__.__name__)


# version
__id__ = "$Id$"

#  End of file 
