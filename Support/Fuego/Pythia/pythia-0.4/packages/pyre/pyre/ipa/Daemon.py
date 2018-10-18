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

from pyre.applications.Daemon import Daemon as BaseDaemon


class Daemon(BaseDaemon):


    def __init__(self, name=None):
        if name is None:
            name = "ipad"
            
        BaseDaemon.__init__(self, name)
        self._userManager = None
        self._sessionManager = None
        return


    def _configure(self):
        BaseDaemon._configure(self)
        self._userManager = self.inventory.userManager
        self._userManager.db = self.inventory.db
        
        self._sessionManager = self.inventory.sessionManager
        self._sessionManager.port = self.inventory.port
        self._sessionManager.userManager = self._userManager
        return


    def _init(self, parent):
        BaseDaemon._init(self, parent)
        return

    
    def _describe(self, pidfile):
        pidfile.write('info = {\n')
        pidfile.write('  "pid": %d,\n' % self._pid)
        pidfile.write('  "port": %d,\n' % self._sessionManager.port)
        pidfile.write('  "passwd": "%s",\n' % self._userManager.db)
        pidfile.write('  }\n')
        return


    def _registerSignalHandlers(self):
        import signal
        signal.signal(signal.SIGHUP, self._userManager.reload)
        return 


    def _serve(self):
        self._sessionManager.serve()
        return


    class Inventory(BaseDaemon.Inventory):

        import pyre.facilities
        import pyre.properties
        from UserManager import UserManager
        from SessionManager import SessionManager

        inventory = [
            pyre.properties.str("db", default="users.db"),
            pyre.properties.int("port", default=50000),

            pyre.facilities.facility("userManager", default=UserManager()),
            pyre.facilities.facility("sessionManager", default=SessionManager())
            ]


# version
__id__ = "$Id$"

#  End of file 
