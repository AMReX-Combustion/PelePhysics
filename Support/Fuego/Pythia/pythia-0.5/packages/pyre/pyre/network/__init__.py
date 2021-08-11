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


# port monitor

from __future__ import absolute_import

_PORT_RANGE_START = 50000
_PORT_RANGE_END = 64 * 1024


def monitor(mode):
    if mode == "tcp":
        from .PortMonitor import PortMonitor

        return PortMonitor()

    if mode == "udp":
        from .UDPMonitor import UDPMonitor

        return UDPMonitor()

    return None


def connection(mode):
    if mode == "tcp":
        from .Client import Client

        return Client()

    if mode == "udp":
        from .UDPConnection import UDPConnection

        return UDPConnection()

    return None


def server():
    from .PortServer import PortServer

    return PortServer()


def client():
    from .Client import Client

    return Client()


def selector():
    from .Selector import Selector

    return Selector()


# version
__id__ = "$Id$"

#  End of file
