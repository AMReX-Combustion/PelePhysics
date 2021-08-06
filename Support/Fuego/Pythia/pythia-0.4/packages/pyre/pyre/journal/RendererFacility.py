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
from pyre.facilities.Facility import Facility


class RendererFacility(Facility):


    def __init__(self, default=None):

        if default is None:
            from .Renderer import Renderer
            default = Renderer()
            
        Facility.__init__(self, name="renderer", default=default)
        return


# version
__id__ = "$Id$"

# End of file 
