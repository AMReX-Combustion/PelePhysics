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
from .Monitor import Monitor


class PortMonitor(Monitor):

    
    def accept(self):
        """encapsulate accept calls on my socket"""
        return self._socket.accept()


    def install(self, port, host='', maxPort=64*1024, queueSize=10):
        """attempt to bind me to the specified port"""

        import socket
        self._install(socket.SOCK_STREAM, port, host, maxPort)
        self._socket.listen(queueSize)

        return self._port, self._socket


# version
__id__ = "$Id$"

# End of file
