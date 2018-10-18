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


from pyre.applications.Application import Application


class RemoteJournalApp(Application):


    def run(self, *args, **kwds):
        mode = self.inventory.mode

        if mode == "server":
            self._server()
        elif mode == "client":
            self._client()
        elif mode == "both":
            self._both()
        else:
            import journal
            firewall = journal.firewall(self.name)
            firewall.log("unknown mode %r" % mode)
            
        return


    def __init__(self):
        Application.__init__(self, "remote")

        self.delay = 0
        self.port = 0
        self.host = "localhost"
        self.method = "udp"

        return


    def _init(self, parent):
        self.port = self.inventory.port
        self.host = self.inventory.host
        self.delay = self.inventory.delay
        return


    def _server(self):
        if not self.port:
            import journal
            journal.error(self.name).log("server: invalid port '%d'" % self.port)
            return
        
        from pyre.journal.Server import Server
        server = Server()
        server.port = self.port
        server.method = self.method
        server.init(self)

        print "server: installed daemon on port '%d'" % server.port
        server.serve()

        return


    def _client(self):
        import journal

        if self.port:

            print "client: sending messages to %r:%d using %r" % (
                self.host, self.port, self.method)
            
            journal.remote(self.port, self.host, self.method)

        for idx in range(5):
            info = journal.debug("test-%02d" % idx).activate()
            info.log("%02d: this is a test" % idx)
        
        return


    def _both(self):
        print "executing both client and server"

        import os

        child = os.fork()
        if child == 0:
            # delay
            import select
            select.select([], [], [], self.delay)
            self._client()
        elif child > 0:
            self._server()
        else:
            import journal
            journal.error(self.name).log("error %d from fork" % child)

        return
            

    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = (
            # host where the journal server is located
            pyre.properties.str("host", default="localhost"),

            # port used by the remote server; 0 means do not use remote facilities
            pyre.properties.int("port", default=0),
            
            # connection mode used to access remote server
            pyre.properties.str(
                "method", default="udp",
                validator=pyre.properties.choice(["udp"])
                ),
            
            # initial delay before launching client
            pyre.properties.int("delay", default=1),
            
            # mode
            pyre.properties.str(
                "mode", default="both",
                validator=pyre.properties.choice(["server", "client", "both"])
                ),
            )


# main

if __name__ == "__main__":
    app = RemoteJournalApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
