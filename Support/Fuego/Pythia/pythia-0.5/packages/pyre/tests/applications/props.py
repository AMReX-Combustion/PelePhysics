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

from pyre.applications.Application import Application


class HelloApp(Application):
    def run(self):
        print("Hello world!")

        self.inventory.bool = True
        self.inventory.str = "hello"
        self.inventory.int = 1
        self.inventory.float = 1.0
        self.inventory.list = [True, "hello", 1, 1.0]
        self.inventory.sequence = "[1-5,8,10]"
        self.inventory.scalar = "1*meter"
        self.inventory.vector = "[1*meter, 1*meter]"

        self.inventory.less = -1.0
        self.inventory.greater = 1.0
        self.inventory.range = 0.9
        self.inventory.choice = "yes"

        print("inventory: %r" % self.state().list())

        return

    def __init__(self):
        super(HelloApp, self).__init__("hello")
        return

    class Inventory(Application.Inventory):

        import pyre.properties
        import pyre.units.SI as si

        inventory = (
            # properties
            pyre.properties.property("property"),
            pyre.properties.bool("bool"),
            pyre.properties.str("str"),
            pyre.properties.int("int"),
            pyre.properties.float("float"),
            pyre.properties.list("list"),
            pyre.properties.sequence("sequence"),
            pyre.properties.dimensional("scalar", default=0.0 * si.meter),
            pyre.properties.dimensional(
                "vector", default=[0.0 * si.meter, 0.0 * si.meter]
            ),
            # validators
            pyre.properties.float("less", validator=pyre.properties.less(0.0)),
            pyre.properties.float(
                "greater", validator=pyre.properties.greater(0.0)
            ),
            pyre.properties.float(
                "range", validator=pyre.properties.range(0.0, 1.0)
            ),
            pyre.properties.str(
                "choice", validator=pyre.properties.choice(["yes", "no"])
            ),
        )


# main
if __name__ == "__main__":
    app = HelloApp()
    app.main()


# version
__id__ = "$Id$"

# End of file
