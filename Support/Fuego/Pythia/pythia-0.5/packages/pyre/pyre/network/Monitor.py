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

from builtins import object

import journal


class Monitor(object):
    def close(self):
        self._socket.close()
        return

    def __init__(self):
        self._port = None
        self._socket = None
        return

    def _install(self, type, port, host="", maxPort=64 * 1024):
        """attempt to bind me to the specified port"""

        if port > maxPort:
            msg = "requested port %d is outside the valid port range" % port
            raise ServerException(msg)

        import socket

        minPort = port
        while port <= maxPort:

            s = socket.socket(socket.AF_INET, type)
            try:
                s.bind((host, port))
                self._port, self._socket = port, s
                self._info.log(
                    "successfully installed at port %d" % self._port
                )
                return

            except socket.error as error:
                number, message = error
                s.close()
                self._info.log(
                    "failed to activate server at port %d: %s"
                    % (port, message)
                )

            port = port + 1

        # no available ports in the range
        msg = "no ports available in the range [%d, %d]" % (minPort, maxPort)
        raise ServerException(msg)

    def _getPort(self):
        return self._port

    def _getSocket(self):
        return self._socket

    # required by select
    def fileno(self):
        return self._socket.fileno()

    port = property(_getPort, None, None, "my port number")
    socket = property(_getSocket, None, None, "my socket object")

    __slots__ = ("_port", "_socket")

    _info = journal.debug("pyre.network.server")


# local exception class


class ServerException(Exception):
    def __init__(self, msg):
        self._message = msg
        return

    def __str__(self):
        return self._message


# version
__id__ = "$Id$"

#  End of file
