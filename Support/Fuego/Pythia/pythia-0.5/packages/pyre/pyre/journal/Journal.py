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

import journal
from pyre.components.Component import Component

from .ChannelFacility import ChannelFacility


class Journal(Component):
    def device(self):
        return self.inventory.device

    def __init__(self):
        Component.__init__(self, "journal", "journal")
        return

    def _init(self, parent):
        import journal

        component = journal.journal()
        device = self.inventory.device.device
        component.device = device

        return

    class Inventory(Component.Inventory):

        from .DeviceFacility import DeviceFacility

        inventory = [
            ChannelFacility(channel)
            for channel in journal.journal().channels()
        ] + [
            DeviceFacility(),
        ]


# version
__id__ = "$Id$"

# End of file
