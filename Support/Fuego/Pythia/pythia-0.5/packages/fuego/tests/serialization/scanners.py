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
def testElementScanner():
    print("Testing element scanner")

    from pyre.chemistry.unpickle.chemkin.scanners.Elements import Elements
    scanner = Elements()

    tests = [
        " ", # Don't send an empty string
        "!",
        "! This is a comment",
        "H2", 
        "H2/1.004/"
        "H2 / 1.004 /"
        ]

    testScanner(scanner, tests)

    return
    

def testSpeciesScanner():
    print("Testing species scanner")

    from pyre.chemistry.unpickle.chemkin.scanners.Species import Species
    scanner = Species()

    tests = [
        " ", # Don't send an empty string
        "!",
        "! This is a comment",
        "H2", 
        "H2O2++",
        "CO2(S)"
        ]

    testScanner(scanner, tests)

    return
    

def testThermoScanner():
    print("Testing thermo scanner")

    from pyre.chemistry.unpickle.chemkin.scanners.Thermo import Thermo
    scanner = Thermo()

    tests = [
        " ", # Don't send an empty string
        "!",
        "! This is a comment",
        "THERMO", 
        "THERMO ALL", 
        ]

    testScanner(scanner, tests)

    return


def testReactionScanner():
    print("Testing reaction scanner")

    from pyre.chemistry.unpickle.chemkin.scanners.Reactions import Reactions
    scanner = Reactions()

    tests = [
        " ", # Don't send an empty string
        "!",
        "! This is a comment",
        "REACTION",
        "CH3+CH3(+M)=C2H6(+M)         9.03E+16 -1.18    653.3      ! 1",
        "Rev / 1.0             2.0                3.0                   /",
        "LOW / 1.0             2.0                3.0                   /",
        "LT / 1.0             2.0   /",
        "RLT / 1.0             2.0   /",
        "SRI / 1.0             2.0   3.0/",
        "SRI / 1.0             2.0   3.0 4.0 5.0 /",
        "TROE / 1.0             2.0   3.0/",
        "TROE / 1.0             2.0   3.0 4.0 /",
        ]

    testScanner(scanner, tests)

    return
    

def testScanner(scanner, tests):
    for test in tests:
        print("    {%s}: %s" % (test, scanner.match(test, 0, 0)))

    return


# main

if __name__ == "__main__":

    testElementScanner()
    testSpeciesScanner()
    testThermoScanner()
    testReactionScanner()


# version
__id__ = "$Id$"

#
# End of file
