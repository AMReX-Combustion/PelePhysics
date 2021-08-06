#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

#factory method

from __future__ import absolute_import
def parser():
    return Parser()


# implementation of the Parser singleton

from pyre.util.Singleton import Singleton


class Parser(Singleton):


    def extend(self, *modules):
        for module in modules:
            self.context.update(module.__dict__)
        return


    def parse(self, string):
        return eval(string, self.context)
    

    def init(self, *args, **kwds):
        self.context = self._initializeContext()
        return


    def _initializeContext(self):
        context = {}

        modules = self._loadModules()
        for module in  modules:
            context.update(module.__dict__)

        return context


    def _loadModules(self):

        from . import SI
        from . import area
        from . import density
        from . import energy
        from . import force
        from . import length
        from . import mass
        from . import power
        from . import pressure
        from . import speed
        from . import substance
        from . import temperature
        from . import time
        from . import volume

        modules = [
            SI,
            area, density, energy, force, length, mass, power, pressure, speed, substance,
            temperature, time, volume
            ]

        return modules
        
    

# version
__id__ = "$Id$"

# End of file 
