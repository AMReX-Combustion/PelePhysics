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


def test():
    pyre.debug.activate("pyre.initialization")

    symbols = [ "He", "Ne", "Ar", "Kr", "Xe", "Rn" ]

    for symbol in symbols:
        element = pyre.chemistry.mechanisms.element(symbol)
        print
        print "%s - atomic weight: %g" % (symbol, element.atomicWeight())
        print "Periodic table: %s" % pyre.handbook.periodicTable().symbol(symbol)

    return


# main

if __name__ == "__main__":
    test()


# version
__id__ = "$Id$"

# End of file
