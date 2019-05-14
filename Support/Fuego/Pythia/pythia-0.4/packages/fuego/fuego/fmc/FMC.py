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

        save          = self.inventory.name
        input         = self.inventory.input
        output        = self.inventory.output
        mechanismFile = self.inventory.mechanism
        thermo        = self.inventory.thermo
        trans         = self.inventory.trans
        #AF
        save_chop   = save.split(".")[0]+"_1.cpp" #self.inventory.name_chop
        save_header = "chemistry_file.H" #save.split(".")[0]+".H" #self.inventory.header

        timer = pyre.monitors.timer("fuego")
        if not input:
            print "Loading '%s'" % (mechanismFile),
        else:
            print "Loading '%s' using '%s' parser" % (mechanismFile, input),

        timer.start()

        mechanism = fuego.serialization.mechanism()
        if thermo:
            mechanism.externalThermoDatabase(thermo)
        if trans:
            mechanism.externalTransDatabase(trans)
        mechanism = fuego.serialization.load(
            filename=mechanismFile, format=input, mechanism=mechanism)
    
        print "... done (%g sec)" % timer.stop()

        timer.reset()
        timer.start()
        print "Converting into '%s' format" % output,
        lines        = fuego.serialization.save(mechanism, output)
        print "... done (%g sec)" % timer.stop()

        print "saving in '%s' (header) and '%s'" % (save_header, save),
        timer.reset()
        timer.start()
        outputFileHeader  = self._openOutput(save_header)
        count_lines = 0
        for line in lines:
            if ('include "chemistry_file.H"') in line:
                line_start_core = count_lines
                break;
            outputFileHeader.write(line)
            outputFileHeader.write('\n')
            count_lines += 1

        outputFile = self._openOutput(save)
        for line in lines[line_start_core:]:
            outputFile.write(line)
            outputFile.write('\n')

        print "... done (%g sec)" % timer.stop()

        return


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
            pyre.properties.str("thermo", default=""),
            pyre.properties.str("trans", default=""),
            pyre.properties.str("input", default=""),
            pyre.properties.str("output", default="c"),
            #pyre.properties.str("output", default="f"),
            ]


# version
__id__ = "$Id$"

#
# End of file
