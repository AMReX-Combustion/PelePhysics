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

        mode = self.inventory.mode
        if mode == "read":
            self.read()
        elif mode == "write":
            registry = self.state()
            self.write(registry)
        else:
            import journal
            journal.firewall(self.name).log("unknown mode '%s'" % mode)

        return


    def read(self):
        import pyre.inventory
        parser = pyre.inventory.parser("pml")

        name = self.inventory.file
        stream = file(name)

        parser.parse(stream)
        registry = parser.document().registry

        self.inventory.file = ""
        self.write(registry)

        return
        
            

    def write(self, registry):
        import pyre.inventory
        renderer = pyre.inventory.renderer()
        renderer.initialize()
        
        document = "\n".join(renderer.render(registry))

        name = self.inventory.file
        if name:
            out = file(name, "w")
        else:
            import sys
            out = sys.stdout

        out.write(document)

        return


    def __init__(self):
        Application.__init__(self, "hello")
        return


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = (
            pyre.properties.str(
                "mode", default="read",
                validator=pyre.properties.choice(["read", "write"])),

            pyre.properties.str("file")
            )


# main
if __name__ == "__main__":
    app = HelloApp()
    app.main()
    

# version
__id__ = "$Id$"

# End of file 
