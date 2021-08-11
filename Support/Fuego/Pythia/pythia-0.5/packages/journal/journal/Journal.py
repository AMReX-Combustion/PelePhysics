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

from builtins import object


class Journal(object):
    def record(self, entry):
        self.device.record(entry)
        return

    def entry(self):
        from .Entry import Entry

        return Entry()

    def channel(self, name, channel=None):
        if channel is None:
            return self._channels.get(name)

        self._channels[name] = channel

        return channel

    def channels(self):
        return list(self._channels.keys())

    def __init__(self, name, device=None):
        self.name = name
        self._channels = {}

        if device is None:
            from .Console import Console

            device = Console()
        self.device = device

        return


# version
__id__ = "$Id$"

#  End of file
