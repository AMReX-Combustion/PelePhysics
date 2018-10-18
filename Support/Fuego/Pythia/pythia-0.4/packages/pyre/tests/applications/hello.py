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


from pyre.applications.Application import Application


class HelloApp(Application):


    def run(self):
        print "Hello %s!" % self.inventory.name
        self._debug.log("Hello world!")
        return
    

    def __init__(self):
        Application.__init__(self, "hello")
        return


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = [
            pyre.properties.str("name", default="world"),
            ]


# main
if __name__ == "__main__":
    app = HelloApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
