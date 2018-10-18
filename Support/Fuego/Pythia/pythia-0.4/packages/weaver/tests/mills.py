#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from pyre.applications.Application import Application


class Stationery(Application):


    def run(self, *args, **kwds):

        filename = self.inventory.name
        language = self.inventory.language

        import weaver
        registry = weaver.registry()
        for language in registry.mills():
            mill = registry.mill(language)
            mill.initialize()

            if filename:
                outfile = file(filename, "w")
            else:
                import sys
                outfile = sys.stdout

            for line in mill.pickle():
                outfile.write(line)
                outfile.write('\n')
        
        return


    def __init__(self):
        Application.__init__(self, "stationery")
        return


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = (
            pyre.properties.str("name"),
            pyre.properties.str("language", default="python"),
            )


# main

if __name__ == "__main__":
    import journal
    # journal.debug("cmdline.parsing").activate()
    # journal.debug("cmdline.configuration").activate()
    # journal.debug("weaver").activate()
    # journal.debug("weaver.mill").activate()
    # journal.debug("weaver.registry").activate()
    # journal.debug("pyre.xml.parsing").activate()

    app = Stationery()
    app.main()


# version
__id__ = "$Id$"

# End of file
