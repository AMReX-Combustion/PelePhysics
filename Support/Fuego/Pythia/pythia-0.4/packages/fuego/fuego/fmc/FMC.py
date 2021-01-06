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
        mech_header = "mechanism.h"

        timer = pyre.monitors.timer("fuego")
        if not input:
            print "\nLoading '%s' as input file" % (mechanismFile)
        else:
            print "\nLoading '%s' using '%s' parser" % (mechanismFile, input)

        timer.start()

        mechanism = fuego.serialization.mechanism()
        if thermo:
            print "Loading '%s' as thermo file" % (thermo)
            mechanism.externalThermoDatabase(thermo)
        if trans:
            print "Loading '%s' as thermo file" % (trans)
            mechanism.externalTransDatabase(trans)
        mechanism = fuego.serialization.load(
            filename=mechanismFile, format=input, mechanism=mechanism)
    
        print "... done (%g sec)" % timer.stop()

        timer.reset()
        timer.start()
        print "\nConverting into '%s' format" % output
        lines        = fuego.serialization.save(mechanism, output)
        print "... done (%g sec)" % timer.stop()

        print "saving in '%s' '%s' (headers) and '%s'" % (mech_header, save_header, save),
        timer.reset()
        timer.start()
        outputFileHeader  = self._openOutput(save_header)
        count_lines = 0
        for line in lines:
            if ('ifndef MECHANISM_CPP') in line:
                line_start_core = count_lines
                break;
            outputFileHeader.write(line)
            outputFileHeader.write('\n')
            count_lines += 1

        outputFile = self._openOutput(save)
        for line in lines[line_start_core:]:
            if ('#ifndef MECHANISM_h') in line:
                line_start_mech_header = count_lines
                break;
            outputFile.write(line)
            outputFile.write('\n')
            count_lines += 1

        MechHeaderFile = self._openOutput(mech_header)
        for line in lines[line_start_mech_header:]:
            MechHeaderFile.write(line)
            MechHeaderFile.write('\n')

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
