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


from pyre.components.Component import Component


class SessionManager(Component):


    def serve(self):
        while not self._done:
            try:
                self._selector.watch(self.inventory.timeout.value)
            except:
                import sys
                type = sys.exc_type
                info = sys.exc_info()[1]
                self._info.log("unhandled '%s': %s" % (type, info))

        return


    def __init__(self, name=None):
        if name is None:
            name = "session-manager"
        Component.__init__(self, name, facility="sessionManager")

        # public data
        self.port = None
        self.timeout = 30 # seconds
        self.userManager = None

        # private data
        self._done = False
        self._tickets = {}
        self._monitor = None
        self._selector = None

        return


    def _init(self, parent):

        import pyre.network
        self._monitor = pyre.network.monitor("tcp")
        self.port, socket = self._monitor.install(self.port)

        self._info.log("session manager installed on port %d" % self.port)

        self._selector = pyre.network.selector()
        self._selector.notifyWhenIdle(self._idle)
        self._selector.notifyOnInterrupt(self._cleanup)
        self._selector.notifyOnReadReady(socket, self._issueTicket)

        return


    def _issueTicket(self, selector, socket):
        import pickle
        from AuthenticationRequest import AuthenticationRequest

        ticketOnce = self.inventory.ticketOnce

        connection, address = socket.accept()

        inbound = connection.makefile("rb")
        outbound = connection.makefile("wb")
        request = pickle.load(inbound)

        username = request.username
        self._info.log("got a request for username '%s'" % (username))

        # check for a ticketed request
        ticket = request.ticket
        if ticket:
            if not self._tickets.has_key(ticket):
                self._info.log("got bad ticket '%s'" % (ticket))
                pickle.dump(0, outbound)
                return True

            if ticket[:len(username)] != username:
                self._info.log("username does not match ticket '%s'" % (ticket))
                pickle.dump(0, outbound)
                return True

            expiration = self._tickets[ticket]

            # one time tickets are invalidated immediately -- the rest just expire
            # this allows the browser's back key to function for a while
            if ticketOnce:
                del self._tickets[ticket]

            import time
            if time.time() - expiration > 0:
                self._info.log("got expired ticket '%s'" % (ticket))
                pickle.dump(0, outbound)
                return True

            newTicket = self._createTicket(username, ticket[len(username):])
            self._info.log("exchanged good ticket '%s' for '%s'" % (ticket, newTicket))

            pickle.dump(newTicket, outbound)
            return True

        username = request.username
        cleartext = request.password
        password = self.userManager.authenticate(username, cleartext)

        if not password:
            pickle.dump(0, outbound)
            return True
        
        ticket = self._createTicket(username, password)
        pickle.dump(ticket, outbound)

        inbound.close()
        outbound.close()
        connection.close()
        
        return True


    def _idle(self, selector):
        import time
        self._info.log("thump")

        currentTime = time.time()
        duration = self.inventory.ticketDuration

        for ticket, timestamp in self._tickets.items():
            if currentTime - timestamp > duration:
                self._info.log("removed stale ticket '%s'" % ticket)
                del self._tickets[ticket]

        return True


    def _createTicket(self, seed, text):
        import time
        import random
        
        t = list(text)
        random.shuffle(t)
        shuffled = "".join(t)

        ticket = seed + shuffled
        self._tickets[ticket] = time.time() + self.inventory.ticketDuration

        return ticket


    def _cleanup(self, selector):
        self._done = True
        return
        
        
    class Inventory(Component.Inventory):

        import pyre.properties
        from pyre.units.time import second

        inventory = [
            pyre.properties.dimensional("timeout", default=10*second),
            pyre.properties.bool("ticketOnce", default=True),
            pyre.properties.float("ticketDuration", default=5),
            ]


# version
__id__ = "$Id$"

#  End of file 
