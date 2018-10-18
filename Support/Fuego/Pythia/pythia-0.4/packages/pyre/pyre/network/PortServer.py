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


from PortMonitor import PortMonitor


class PortServer(PortMonitor):


    def deactivate(self):
        self._selector.state = False
        return

    
    def activate(self, onConnection, onTimeout, timeout):
        """enter an infinite loop, calling the specified callbacks when appropriate"""

        import pyre.network
        selector = pyre.network.selector

        self._connectionHandler = onConnection
        selector.notifyOnTimeout(self._socket, onTimeout)
        selector.notifyOnReadReady(self._socket, self._onConnection)
        
        self._selector = selector
        selector.watch(timeout)

        return


    def __init__(self):
        PortMonitor.__init__(self)
        self._selector = None
        self._connectionHandler = None
        return


    def _onConnection(self, *args):
        descriptor, address = self.accept()
        self._info.log("connection attempt from '%s' on port %d" % (address, self._port))

        # the callback determines whether we exit or wait some more
        return self._handler(descriptor)



# version
__id__ = "$Id$"

# End of file
