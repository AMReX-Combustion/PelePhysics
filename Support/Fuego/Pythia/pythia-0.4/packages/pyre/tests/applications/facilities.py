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


import pyre.properties
from pyre.components.Component import Component
from pyre.applications.Application import Application


def facility():
    component = TestComponent()
    component.name = "component1"
    component.inventory.property = "script"
    return component


class TestComponent(Component):


    def value(self):
        return self.inventory.property


    def __init__(self):
        Component.__init__(self, "component", "facility")
        return


    class Inventory(Component.Inventory):

        inventory = (
            pyre.properties.str("property", default="value"),
            )


class HelloApp(Application):


    def run(self):
        print "property='%s'" % (self.inventory.facility.value())
        self._debug.log("journal: '%s'" % self.inventory.journal)
        return


    def __init__(self):
        super(HelloApp, self).__init__("hello")
        return


    class Inventory(Application.Inventory):

        inventory = (
            pyre.facilities.facility("facility", default=TestComponent()),
            )


# driver

def test():
    app = HelloApp()
    app.main()
    return


# main
if __name__ == "__main__":
    test()


# version
__id__ = "$Id$"

# End of file 
