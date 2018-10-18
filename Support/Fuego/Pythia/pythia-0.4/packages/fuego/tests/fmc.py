#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.applications.Application import Application


class FMC(Application):


    def run(self):

        import fuego
        import pyre.monitors

        save = self.inventory.name
        input = self.inventory.input
        output = self.inventory.output
        mechanismFile = self.inventory.mechanism

        timer = pyre.monitors.timer("fuego")
        if not input:
            print "Loading '%s'" % (mechanismFile),
        else:
            print "Loading '%s' using '%s' parser" % (mechanismFile, input),

        timer.start()
        mechanism = fuego.serialization.load(mechanismFile, input)
        print "... done (%g sec)" % timer.stop()

        timer.reset()
        timer.start()
        print "Converting into '%s' format" % output,
        lines = fuego.serialization.save(mechanism, output)
        print "... done (%g sec)" % timer.stop()

        print "saving in '%s' ..." % save
        outputFile = self._openOutput(save)
        for line in lines:
            outputFile.write(line)
            outputFile.write('\n')

        return mechanism


    def __init__(self):
        Application.__init__(self, "fmc")
        return


    def _openOutput(self, name):
        if name == "stdout":
            import sys
            return sys.stdout

        return file(name, "w")


    class Inventory(Application.Inventory):

        import pyre.properties

        inventory = [
            pyre.properties.str("name", default="stdout"),
            pyre.properties.str("mechanism", default="GRIMech-3.0.ck2"),
            pyre.properties.str("input", default=""),
            pyre.properties.str("output", default="c"),
            ]


# main

if __name__ == "__main__":
    import journal
    journal.info("fuego").activate()
    journal.debug("fuego").activate()

    app = FMC()
    app.main()

    

# version
__id__ = "$Id$"

#
# End of file
