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


from builtins import object


class UDPConnection(object):
    def sendto(self, message, port, host=""):
        server = (host, port)
        self._socket.sendto(message, server)
        return

    def __init__(self):
        import socket

        localinfo = socket.gethostbyaddr(socket.gethostname())

        self._localhost = localinfo[0]
        self._localip = localinfo[2][0]
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

        return


# version
__id__ = "$Id$"

# End of file
