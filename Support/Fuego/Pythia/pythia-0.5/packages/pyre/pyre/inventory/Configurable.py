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


from __future__ import absolute_import

from builtins import object


class Configurable(object):
    def defaults(self):
        self.inventory.defaults()
        self._defaults()
        return

    def configure(self, registry):
        report = self.inventory.configure(registry)
        self._configure()
        return report

    def init(self, parent=None):
        self.preinit()
        self.inventory.init(self)
        self._init(parent)
        self.postinit()
        return

    def fini(self):
        self._fini()
        self.inventory.fini()
        return

    def preinit(self):
        return

    def postinit(self):
        return

    def state(self, registry=None):
        if registry is None:
            from .Registry import Registry

            registry = Registry(self.name)

        return self.inventory.state(registry)

    def __init__(self, name):
        self.name = name
        self.inventory = self.Inventory()
        return

    def _defaults(self):
        return

    def _init(self, parent):
        return

    def _configure(self):
        return

    def _fini(self):
        return

    # inventory
    from .Inventory import Inventory


# version
__id__ = "$Id$"

# End of file
