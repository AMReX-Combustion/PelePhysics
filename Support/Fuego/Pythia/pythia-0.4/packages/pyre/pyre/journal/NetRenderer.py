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

from pyre.components.Component import Component


class NetRenderer(Component):


    def __init__(self, name="renderer"):
        Component.__init__(self, name, "net-renderer")
        self.renderer = None
        return


    def _init(self, parent):
        from journal.NetRenderer import NetRenderer
        renderer = NetRenderer()
        self.renderer = renderer
        
        return renderer



# version
__id__ = "$Id$"

# End of file 
