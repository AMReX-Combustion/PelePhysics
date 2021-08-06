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

from __future__ import print_function
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
        mech_header = "mechanism.H"

        timer = pyre.monitors.timer("fuego")
        if not input:
            print("\nLoading '%s' as input file" % (mechanismFile))
        else:
            print("\nLoading '%s' using '%s' parser" % (mechanismFile, input))

        timer.start()

        mechanism = fuego.serialization.mechanism()
        if thermo:
            print("Loading '%s' as thermo file" % (thermo))
            mechanism.externalThermoDatabase(thermo)
        if trans:
            print("Loading '%s' as thermo file" % (trans))
            mechanism.externalTransDatabase(trans)
        mechanism = fuego.serialization.load(
            filename=mechanismFile, format=input, mechanism=mechanism)
    
        print("... done (%g sec)" % timer.stop())

        timer.reset()
        timer.start()
        print("\nConverting into '%s' format" % output)
        lines        = fuego.serialization.save(mechanism, output)
        print("... done (%g sec)" % timer.stop())

        print("saving in '%s' (header) and '%s'" % (mech_header, save), end=' ')
        timer.reset()
        timer.start()
        outputFile = self._openOutput(save)
        count_lines = 0
        for line in lines:
            if ('ifndef MECHANISM_H') in line:
                line_start_core = count_lines
                break;
            outputFile.write(line)
            outputFile.write('\n')
            count_lines += 1

        outputFileHeader  = self._openOutput(mech_header)
        for line in lines[line_start_core:]:
            outputFileHeader.write(line)
            outputFileHeader.write('\n')

        print("... done (%g sec)" % timer.stop())

        return


    def __init__(self):
        Application.__init__(self, "fmc")
        return


    def _openOutput(self, name):
        if name == "stdout":
            import sys
            return sys.stdout

        return open(name, "w")


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
