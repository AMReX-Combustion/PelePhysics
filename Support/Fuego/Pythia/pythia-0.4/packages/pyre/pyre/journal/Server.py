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

from Journal import Journal


class Server(Journal):


    def serve(self, timeout=10):
        import pickle
        import select
        from journal.Entry import Entry

        import journal
        theJournal = journal.journal()
        file = self.monitor.socket.makefile("rb")

        while 1:
            reads, writes, excepts = select.select([file], [], [], timeout)
            if not reads:
                self._idle()
                continue

            entry = pickle.load(file)
            theJournal.record(entry)

        return


    def __init__(self):
        Journal.__init__(self)
        self.port = 50000
        self.method = "udp"
        self.monitor = None
        return


    def _init(self, parent):
        Journal._init(self, parent)
        self.monitor = self._setupMonitor(self.port, self.method)
        return


    def _setupMonitor(self, port, method):
        import pyre.network
        monitor = pyre.network.monitor(method)
        monitor.install(port)

        return monitor


    def _idle(self):
        return


# version
__id__ = "$Id$"

#  End of file 
