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

from __future__ import print_function
from __future__ import absolute_import
from pyre.components.Component import Component


class Application(Component):


    def main(self, *args, **kwds):
        self.defaults()

        registry = self.prime()
        self.configure(registry)

        self.init()

        self.execute(*args, **kwds)

        self.fini()

        return


    def prime(self):
        parser = self.inventory.commandlineParser
        self._registry.orphans += parser.parse(root=self._registry)
        return self._registry


    def execute(self, *args, **kwds):
        if self._registry.help:
            self.help()
        elif self._registry.unknown:
            self.usage()
        else:
            self.run(*args, **kwds)

        return


    def help(self):
        return


    def usage(self):
        print("unkown arguments:", self._registry.unknown)
        return


    def __init__(self, name):
        Component.__init__(self, name, "application")

        import sys
        self.filename = sys.argv[0]
        self._registry = self._createRegistry()

        return


    def _createRegistry(self):
        from .ConfigurationRegistry import ConfigurationRegistry
        return ConfigurationRegistry(self.name)


    class Inventory(Component.Inventory):

        import pyre.journal
        import pyre.facilities
        from .CommandlineParser import CommandlineParser

        inventory = [
            pyre.journal.journal(),
            pyre.facilities.facility("commandlineParser", default=CommandlineParser())
            ]


# version
__id__ = "$Id$"

# End of file 
