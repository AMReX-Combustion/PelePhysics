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

from .Indenter import Indenter
from .Stationery import Stationery


class Mill(Stationery, Indenter):
    def pickle(self, document=None):

        self._begin()
        if document:
            self._renderDocument(document)
        self._end()

        return self._rep

    def initialize(self, options=None):
        import weaver.config.unpickle

        self._debug.line("reading user configuration")
        userOptions = weaver.config.unpickle.readUserConfiguration()

        if userOptions:

            self._debug.line("user options: %s" % repr(userOptions))
            # pull out the defaults
            self._debug.line("extracting default options")
            defaults = userOptions.facilities.get("default")
            if defaults:
                self._debug.line("found default options")
                self.inventory.configure(defaults)

            # pull out the mill specific section
            self._debug.line("extracting mill-specific options")
            mill = userOptions.facilities.get(self.names[0])
            if mill:
                self._debug.line("found mill-specific options")
                self.inventory.configure(mill)

        # override with the ones from the argument list
        if options:
            self.inventory.configure(options)

        self._debug.log()

        return

    def __init__(self):
        Stationery.__init__(self, "mill")
        Indenter.__init__(self)

        self._rep = []

        return

    def _begin(self):
        # self._rep = self.header()
        return

    def _end(self):
        self._versionId()
        self._rep += ["", self.line(self._timestamp())]
        self._rep += self.footer()
        return

    def _separator(self):
        self._rep.append(self.line(self.separator()))
        return

    def _versionId(self):
        format = self.inventory.versionId

        if format:
            self._rep += ["", self.line(" version"), self.line(format)]

        return

    def _timestamp(self):
        format = self.inventory.timestampLine

        if format:
            import time

            timestamp = format % (self.__class__.__name__, time.asctime())
            return timestamp

        return ""

    def _write(self, text=""):
        if text:
            self._rep.append(self._margin + text)
        else:
            self._rep.append("")
        return

    # properties
    class Inventory(Stationery.Inventory):

        import pyre.properties

        inventory = (
            pyre.properties.property(
                "timestampLine", default=" Generated automatically by %s on %s"
            ),
            pyre.properties.property("versionId", default=" $" + "Id" + "$"),
        )


# version
__id__ = "$Id$"

#  End of file
