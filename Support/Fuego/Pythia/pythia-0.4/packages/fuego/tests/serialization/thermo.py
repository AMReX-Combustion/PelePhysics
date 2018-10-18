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


def load(options):

    pyre.support.debug.activate("pyre.chemistry.serialization")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-scanner")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-parser")
    # pyre.support.debug.activate("pyre.chemistry.chemkin-tokenizer")
    
    file = options["--file"]
    format = options["--format"]

    timer = pyre.timers.timer("chemkin-unpickle")
    print "Loading '%s'; assuming '%s' format" % (file, format)

    timer.start()
    db = pyre.chemistry.serialization.loadThermoDatabase(file, format)
    timer.stop()

    print "Done loading '%s' (%g sec)" % (file, timer.read())

    if db:
        db.dump()

    return db


# usage

def usage(program):

    formats = pyre.chemistry.serialization.unpicklers().registered()
    supported = "|".join(formats)
    
    print "Usage: %s [options ...]" % program
    print "Options: (default values in brackets)"
    print "    --file=<mechanism filename> [%s]" % defaults["--file"]
    print "    --format=<%s> [%s]" % (supported, defaults["--format"])
    print
    return
        

defaults = {
    "--file": "therm.dat",
    "--format": "chemkin"
    }

# main

if __name__ == "__main__":

    options = pyre.applications.main(defaults, usage)
    db = load(options)
    

# version
__id__ = "$Id$"

#
# End of file
