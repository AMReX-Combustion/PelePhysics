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

from .Device import Device


class Remote(Device):
    def createDevice(self):

        host = self.inventory.host
        port = self.inventory.port
        mode = self.inventory.mode

        from . import journal

        return journal.remote(port=port, host=host, mode=mode)

    def __init__(self):
        Device.__init__(self, "remote")
        return

    class Inventory(Device.Inventory):

        import pyre.properties

        from .NetRenderer import NetRenderer
        from .RendererFacility import RendererFacility

        inventory = (
            pyre.properties.str("host", default="localhost"),
            pyre.properties.int(
                "port",
                default=50000,
                validator=pyre.properties.range(1024 + 1, 64 * 1024 - 1),
            ),
            pyre.properties.str(
                "mode",
                default="udp",
                validator=pyre.properties.choice(["udp"]),
            ),
            RendererFacility(default=NetRenderer()),
        )


# version
__id__ = "$Id$"

# End of file
