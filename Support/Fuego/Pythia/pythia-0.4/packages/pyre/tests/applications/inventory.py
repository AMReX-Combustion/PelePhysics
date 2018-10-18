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

        import pyre.properties

        inventory = (
            pyre.properties.str("property", default="value"),
            )


class HelloApp(Application):


    def run(self):

        self.listProperties()
        self.listComponents()

        print "facility.property = '%s'" % self.inventory.facility.value()

        return


    def preinit(self):
        self.inventory.name = "Michael"
        self.Inventory.name.default = "Keri"
        return


    def listProperties(self):
        print "properties: {"
        for name, prop in self.inventory.properties().iteritems():
            default = prop.default
            tp = prop.type
            value = self.inventory.__getattribute__(name)
            print "  (name='%s', type='%s', default='%s', value='%s')" % (name, tp, default, value)
        print "  }"
        return


    def listComponents(self):
        print "components: {"
        for name, facility in self.inventory.facilities().iteritems():
            default = facility.default
            value = self.inventory.__getattribute__(name)
            print "  (facility='%s', default='%s', component='%s')" % (
                name, default.name, value.name)
        print "  }"
        return


    def __init__(self):
        Application.__init__(self, "hello")
        return


    class Inventory(Application.Inventory):

        import pyre.properties
        import pyre.facilities

        inventory = [
            pyre.properties.str("name", default="world"),
            pyre.facilities.facility("facility", default=TestComponent())
            ]


# main
if __name__ == "__main__":
    app = HelloApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
