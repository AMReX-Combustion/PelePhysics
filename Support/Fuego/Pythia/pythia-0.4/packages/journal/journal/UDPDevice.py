#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


import pickle
from Device import Device


class UDPDevice(Device):


    def record(self, entry):
        str = self.renderer.render(entry)
        self._socket.sendto(pickle.dumps(str), self._server)

        return


    def __init__(self, port, host=''):
        import socket
        from NetRenderer import NetRenderer

        Device.__init__(self, NetRenderer())

        self._server = (host, port)
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

        return


# version
__id__ = "$Id$"

# End of file
