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

import traceback
from builtins import object

import journal

from .Entry import Entry


class Diagnostic(object):
    def line(self, message):
        if not self._state:
            return

        self._entry.line(message)
        return self

    def log(self, message=None):
        if not self._state:
            return

        if message is not None:
            self._entry.line(message)

        stackDepth = -2
        stackTrace = traceback.extract_stack()
        file, line, function, src = stackTrace[stackDepth]

        meta = self._entry.meta
        meta["category"] = self.name
        meta["facility"] = self.facility
        meta["filename"] = file
        meta["function"] = function
        meta["line"] = line
        meta["src"] = src
        meta["stack-trace"] = stackTrace[: stackDepth + 1]

        journal.journal().record(self._entry)

        if self._fatal:
            raise self.Fatal

        self._entry = Entry()
        return self

    def activate(self):
        self._state = True
        return self

    def deactivate(self):
        self._state = False
        return self

    def flip(self):
        self._state ^= True
        return self

    def __init__(self, name, facility, defaultState, fatal=False):
        self.name = name
        self.facility = facility

        self._entry = Entry()
        self._state = defaultState
        self._fatal = fatal

        return

    def _getState(self):
        return self._state

    def _setState(self, state):
        self._state = state
        return

    state = property(_getState, _setState, None, "")

    class Fatal(Exception):
        def __str__(self):
            return "fatal diagnostic"


# version
__id__ = "$Id$"

#  End of file
