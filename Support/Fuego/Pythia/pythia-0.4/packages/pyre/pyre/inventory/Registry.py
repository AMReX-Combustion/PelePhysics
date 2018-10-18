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


class Registry(object):


    def open(self, name):
        self._debug.line("%s: looking for node '%s'" % (self.name, name))

        try:
            node = self.facilities[name]
            self._debug.line("%s: found node '%s'" % (self.name, name))
        except KeyError:
            node = Registry(name)
            self.facilities[name] = node
            self._debug.line("%s: created node '%s'" % (self.name, name))

        self._debug.log()
        return node


    def attach(self, node):
        self.facilities[node.name] = node
        return


    def set(self, name, value):
        self.properties[name] = value
        return


    def list(self):

        listing = [
            ("%s.%s" % (self.name, name), value)
            for name, value in self.properties.iteritems()
            ]

        listing += [
            ("%s.%s" % (self.name, name), value)
            for facility in self.facilities.itervalues()
            for name, value in facility.list() 
            ]

        return listing


    def __init__(self, name):
        self.name = name
        self.properties = {}
        self.facilities = {}

        import journal
        self._debug = journal.debug("registry")
        
        return

# version
__id__ = "$Id$"

# End of file 
