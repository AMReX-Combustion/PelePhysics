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
from pyre.applications.Toggle import Toggle
from pyre.applications.Resource import Resource


class Timer(Resource, Toggle):


    def start(self):
        if self.state:
            self._start = time.clock()

        return self


    def stop(self):
        now = time.clock()
        self._accumulatedTime += now - self._start
        self._start = now

        return self._accumulatedTime


    def lap(self):
        if self.state:
            now = time.clock()
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

        self._start = time.clock()
        self._accumulatedTime = 0

        return


# version
__id__ = "$Id$"

#  End of file 
