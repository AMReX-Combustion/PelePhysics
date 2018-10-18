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


class NetRenderer(object):


    def render(self, entry):
        entry.meta["ip"] = self._localip
        entry.meta["host"] = self._localhost
        return entry


    def __init__(self):
        import socket

        localinfo = socket.gethostbyaddr(socket.gethostname())

        self._localhost = localinfo[0]
        self._localip = localinfo[2][0]

        return
            

# version
__id__ = "$Id$"

# End of file 
