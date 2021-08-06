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
import pickle
import journal
import pyre.network
from pyre.inventory.Configurable import Configurable


class Session(Configurable):


    def login(self, username, cleartext):
        self._info.log("login request for user '%s'" % username)

        from .AuthenticationRequest import AuthenticationRequest
        request = AuthenticationRequest(username, password=cleartext)
        ticket = self._authenticate(request)
        return ticket


    def resume(self, username, ticket):
        self._info.log("ticketed request for user '%s':'%s'" % (username, ticket))

        from .AuthenticationRequest import AuthenticationRequest
        request = AuthenticationRequest(username, ticket=ticket)
        ticket = self._authenticate(request)
        return ticket


    def logout(self, username, ticket):
        return ""


    def info(self):
        self._info.log("host=%s, port=%d" % (self.inventory.host, self.inventory.port))
        return
        

    def __init__(self, name=None):
        if name is None:
            name = "session"
        Configurable.__init__(self, name)
        self._info = journal.info(name)
        return


    def _authenticate(self, request):

        host = self.inventory.host
        port = self.inventory.port
        method = self.inventory.method

        client = pyre.network.connection(method)

        try:
            client.connect(host, port)
        except client.ClientException:
            self._info.log("unable to connect to authentication server")
            return ""

        client.send(request.encode())

        inbound = client.socket().makefile("rb")
        ticket = pickle.load(inbound)
        self._info.log("received ticket '%s'" % ticket)

        return ticket


    # properties
    class Inventory(Configurable.Inventory):


        import pyre.properties


        inventory = (
            pyre.properties.str('host', 'localhost'),
            pyre.properties.int('port', 50000),
            pyre.properties.str('method', 'tcp'),
            pyre.properties.str('mailto', ''),
            )



# version
__id__ = "$Id$"

#  End of file 
