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
from pyre.components.Component import Component


class Renderer(Component):


    def __init__(self, name="renderer"):
        Component.__init__(self, name, "renderer")
        self.renderer = None
        return


    def _init(self, parent):
        from journal.Renderer import Renderer
        renderer = Renderer()

        renderer.header = self.inventory.header
        renderer.footer = self.inventory.footer
        renderer.format = self.inventory.format

        self.renderer = renderer
        
        return renderer


    class Inventory(Component.Inventory):


        import pyre.properties

        inventory = (
            pyre.properties.str(
                "header",
                default=" >> %(category)s(%(facility)s) -- %(filename)s:%(line)s:%(function)s"),
            pyre.properties.str("footer", default=""),
            pyre.properties.str("format", default=" -- %s"),
            )


# version
__id__ = "$Id$"

# End of file 
