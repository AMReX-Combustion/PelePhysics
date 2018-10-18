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

import journal


class Client(object):


    def socket(self):
        return self._socket


    def connect(self, host, port):
        import socket
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        self._info.log("attempting connection to '%s', port '%s'" % (host, port))
        try:
            self._socket.connect((host, port))

        except socket.error, error:
            number, reason = error
            raise self.ClientException(host, port, number, reason)

        self._info.log("successful connection to '%s', port %d" % (host, port))

        return


    def close(self):
        return self._socket.close()


    def send(self, bytes):
        return self._socket.send(bytes)


    def recv(self, maxBytes):
        return self._socket.recv(maxBytes)


    def __init__(self):
        self._socket = None
        return


    # static members
    _info = journal.debug("pyre.network.client")


    # local exception class
    
    class ClientException(Exception):

        def __init__(self, host, port, errno, reason):
            self.host = host
            self.port = port
            self.errno = errno
            self.reason = reason
            return

        def __str__(self):
            msg = "error %d connecting to '%s', port %d: %s" % (
                self.errno, self.host, self.port, self.reason)
            return msg


# version
__id__ = "$Id$"

#  End of file 
