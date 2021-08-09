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
from . import journal
from pyre.applications.Daemon import Daemon as BaseDaemon


class Daemon(BaseDaemon):


    def __init__(self, name="journald"):
        BaseDaemon.__init__(self, name)
        self._server = None
        return


    def _init(self, parent):
        BaseDaemon._init(self, parent)

        self._server = self.inventory.journal
        self._logfile = self.inventory.logfile

        return


    def _describe(self, pidfile):
        pidfile.write('info = {\n')
        pidfile.write('  "pid": %d,\n' % self._pid)
        pidfile.write('  "port": %d,\n' % self._server.port)
        pidfile.write('  "logfile": "%s",\n' % self._logfile)
        pidfile.write('  }\n')
        return


    def _serve(self):
        self._server.serve(self.inventory.timeout)
        return


    class Inventory(BaseDaemon.Inventory):


        import pyre.properties
        from .Server import Server


        inventory = (
            pyre.properties.int("port", default=50000),
            pyre.properties.int("timeout", default=30),

            pyre.journal.journal(Server()),
            )


# version
__id__ = "$Id$"

#  End of file 
