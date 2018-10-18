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

import os
from Application import Application


class Daemon(Application):


    def run(self, *args, **kwds):
        self._debug.log("launching daemon...")

        if self.inventory.debug:
            print "debug mode..."
            import os
            self._pid = os.getpid()
            self.serve()
            return True

        import pyre.util
        return pyre.util.spawn(self.done, self.respawn)


    def done(self, pid):
        return True


    def respawn(self, pid):
        self._debug.log("respawning daemon...")

        os.chdir("/")
        os.setsid()
        os.umask(0)

        import pyre.util
        pyre.util.spawn(self.exit, self.daemon)
        return


    def exit(self, pid):
        import sys
        sys.exit(0)


    def daemon(self, pid):
        self._debug.log("daemon running...")

        home = self.inventory.home
        os.close(2)
        os.close(1)
        os.close(0)

        os.chdir(home)

        self._pid = pid
        self.serve()
        
        return
    

    def serve(self):
        self._debug.log("registering signal handlers")
        self._registerSignalHandlers()

        self._debug.log("saving process information")
        self._saveInfo()

        self._debug.log("entring main daemon loop")
        self._serve()
        return


    def __init__(self, name):
        Application.__init__(self, name)
        self._pid = None
        return


    def _defaults(self):
        import pyre.util
        tmpdir = pyre.util.tmp()
        self.Inventory.logfile = os.path.join(tmpdir, "%s.log" % self.name)
        self.Inventory.pidfile = os.path.join(tmpdir, "%s.pid" % self.name)
        return


    def _configure(self):
        import pyre.journal
        device = pyre.journal.file()
        device.inventory.name = self.inventory.logfile
        self.inventory.journal.inventory.device = device
        return


    def _saveInfo(self):
        try:
            filename = self.inventory.pidfile
            pidfile = file(filename, "w")
        except:
            filename = self.Inventory.pidfile
            pidfile = file(filename, "w")

        self._describe(pidfile)
        pidfile.close()

        return


    def _describe(self, pidfile):
        return


    def _registerSignalHandlers(self):
        return


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = [
            pyre.properties.bool("debug", default=False),
            pyre.properties.str("home", default="/"),
            pyre.properties.str("pidfile", default="/tmp/daemon.pid"),
            pyre.properties.str("logfile", default="/tmp/daemon.log"),
            ]


# version
__id__ = "$Id$"

#  End of file 
