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

from __future__ import print_function
import pyre


def load(options):

    pyre.support.debug.activate("pyre.chemistry.serialization")
    pyre.support.debug.activate("fuego.serialization")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-scanner")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-parser")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-tokenizer")
    
    file = options["--file"]
    format = options["--format"]

    timer = pyre.timers.timer("chemkin-unpickle")
    if not format:
        print("Loading '%s'" % (file))
    else:
        print("Loading '%s'; assuming '%s' format" % (file, format))

    timer.start()
    mechanism = pyre.chemistry.serialization.load(file, format)
    timer.stop()

    print("Done loading '%s' (%g sec)" % (file, timer.read()))

    if mechanism:
        mechanism.dump()

    return mechanism


# usage

def usage(program):

    formats = pyre.chemistry.serialization.unpicklers().registered()
    supported = "|".join(formats)
    
    print("Usage: %s [options ...]" % program)
    print("Options: (default values in brackets)")
    print("    --file=<mechanism filename> [%s]" % defaults["--file"])
    print("    --format=<%s> [%s]" % (supported, defaults["--format"]))
    print()
    return
        

defaults = {
    #"--file": "benzene.ck2",
    #"--file": "HMX-CIT-20010831.ck2",
    "--file": "GRIMech-3.0.ck2",
    "--format": ""
    }

# main

if __name__ == "__main__":

    options = pyre.applications.main(defaults, usage)
    mechanism = load(options)
    

# version
__id__ = "$Id$"

#
# End of file
