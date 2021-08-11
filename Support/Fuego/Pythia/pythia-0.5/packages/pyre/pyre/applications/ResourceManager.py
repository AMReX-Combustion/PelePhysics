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


from __future__ import print_function

from builtins import object


class ResourceManager(object):
    def find(self, name):
        return self._registry.get(name)

    def manage(self, resource, id, aliases=[]):
        existing = self._resources.get(resource)
        if existing:
            print(" *++* existing resource:", repr(existing))
            pyre.journal.warning()

        self._resources[resource] = tuple([id] + aliases)
        self._registry[id] = resource
        for alias in aliases:
            self._registry[alias] = resource

        return

    def name(self):
        return self._name

    def resources(self):
        return list(self._resources.keys())

    def __init__(self, name):
        self._name = name

        self._registry = {}
        self._resources = {}

        return


# version
__id__ = "$Id: ResourceManager.py,v 1.1 2003/06/02 19:27:40 aivazis Exp $"

#  End of file
