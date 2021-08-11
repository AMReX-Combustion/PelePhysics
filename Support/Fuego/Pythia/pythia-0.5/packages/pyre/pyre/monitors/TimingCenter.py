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

# factory

from __future__ import absolute_import

from pyre.applications.ResourceManager import ResourceManager


def timingCenter():
    global _theTimingCenter
    if _theTimingCenter is None:
        _theTimingCenter = TimingCenter()

    return _theTimingCenter


# implementation


class TimingCenter(ResourceManager):
    def timer(self, name):
        timer = self.find(name)
        if not timer:
            from .Timer import Timer

            timer = Timer(name).activate()
            self.manage(timer, name)

        return timer

    def __init__(self):
        ResourceManager.__init__(self, "timers")
        return


# the instance

_theTimingCenter = None


# version
__id__ = "$Id$"

#  End of file
