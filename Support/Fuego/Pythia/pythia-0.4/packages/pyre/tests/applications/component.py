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
        print("Hello %s!" % self.inventory.person)

        registry = self.state()

        import pyre.inventory
        renderer = pyre.inventory.renderer()
        renderer.initialize()
        document = "\n".join(renderer.render(registry))
        print(document)

        return


    def __init__(self):
        Application.__init__(self, "hello")
        return


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = (
            pyre.properties.str("person", default="world"),
            )


# main
if __name__ == "__main__":
    app = HelloApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
