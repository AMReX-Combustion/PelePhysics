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

from Monitor import Monitor


class UDPMonitor(Monitor):

    
    def install(self, port, host='', maxPort=64*1024):
        """attempt to bind me to the specified port"""

        import socket
        self._install(socket.SOCK_DGRAM, port, host, maxPort)

        return self._port, self._socket


# version
__id__ = "$Id$"

#  End of file 
