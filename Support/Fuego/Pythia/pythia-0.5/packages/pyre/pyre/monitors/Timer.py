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


import time

from pyre.applications.Resource import Resource
from pyre.applications.Toggle import Toggle


class Timer(Resource, Toggle):
    def start(self):
        if self.state:
            self._start = self.get_clock()

        return self

    def stop(self):
        now = self.get_clock()
        self._accumulatedTime += now - self._start
        self._start = now

        return self._accumulatedTime

    def lap(self):
        if self.state:
            now = self.get_clock()
            return self._accumulatedTime + (now - self._start)

        return 0

    def read(self):
        return self._accumulatedTime

    def reset(self):
        self._accumulatedTime = 0
        return self

    def __init__(self, name):
        Resource.__init__(self, name)
        Toggle.__init__(self)

        self._start = self.get_clock()
        self._accumulatedTime = 0

        return

    def get_clock(self):
        try:
            return time.clock()
        except AttributeError:
            return time.process_time()


# version
__id__ = "$Id$"

#  End of file
