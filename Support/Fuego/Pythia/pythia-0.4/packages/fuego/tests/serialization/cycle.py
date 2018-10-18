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


def load(filename):
    mechanism = pyre.chemistry.serialization.load(filename)
    return mechanism


def save(filename, mechanism):

    file = open(filename, "w")
    lines = pyre.chemistry.serialization.save(mechanism)

    for line in lines:
        file.write(line+'\n')

    file.close()

    return


def cycle(options):

    file = options["--file"]
    format = options["--format"]

    take_1 = load(file)
    # take_1.dump()
    save("take_1.ck2", take_1)

    take_2 = load("take_1.ck2")
    # take_2.dump()
    save_2 ="take_2.ck2"
    save("take_2.ck2", take_2)

    return take_2


# usage

def usage(program):
    print "Usage: %s [options ...]" % program
    print "Options: (default values in brackets)"
    print "    --file=<mechanism filename> [%s]" % defaults["--file"]
    print "    --format=<chemkin|ckml> [%s]" % defaults["--file"]
    print
    return
        

defaults = {
    "--file": "GRIMech-3.0.ck2",
    "--format": "chemkin"
    }

# main

if __name__ == "__main__":

    import pyre

    options = pyre.applications.main(defaults, usage)
    mechanism = cycle(options)
    

# version
__id__ = "$Id$"

#
# End of file
