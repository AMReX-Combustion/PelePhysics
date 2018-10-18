#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        (C) 1998-2001  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

import pyre


def save(options):

    save = options["--file"]
    input = options["--input"]
    output = options["--output"]
    mechanismFile = options["--mechanism"]

    timer = pyre.timers.timer("chemkin-unpickle")
    if not input:
        print "Loading '%s'" % (mechanismFile),
    else:
        print "Loading '%s' using '%s' parser" % (mechanismFile, input),

    timer.start()
    mechanism = pyre.chemistry.serialization.load(mechanismFile, input)
    print "... done (%g sec)" % timer.stop()

    timer.reset()
    timer.start()
    print "Converting into '%s' format" % output,
    lines = pyre.chemistry.serialization.save(mechanism, output)
    print "... done (%g sec)" % timer.stop()

    print "saving in '%s' ..." % save
    outputFile = openOutput(save)
    for line in lines:
        outputFile.write(line)
        outputFile.write('\n')
        
    return mechanism


def openOutput(file):
    if file == "stdout":
        import sys
        return sys.stdout
    return open(file, "w")


# usage

def usage(program):

    inputs = '|'.join(pyre.chemistry.serialization.unpicklers().registered())
    outputs = '|'.join(pyre.chemistry.serialization.picklers().registered())
    
    print "Usage: %s [options ...]" % program
    print "Options: (default values in brackets)"
    print "    --file=<output filename> [%s]" % defaults["--file"]
    print "    --mechanism=<mechanism filename> [%s]" % defaults["--mechanism"]
    print "    --input=<%s> [%s]" % (inputs, defaults["--input"])
    print "    --output=<%s> [%s]" % (outputs, defaults["--output"])
    print
    return
        

defaults = {
    "--file": "stdout",
    "--mechanism": "GRIMech-3.0.ck2",
    "--input": "",
    "--output": "c"
    }

# main

if __name__ == "__main__":

    options = pyre.applications.main(defaults, usage)
    mechanism = save(options)
    

# version
__id__ = "$Id$"

#
# End of file
