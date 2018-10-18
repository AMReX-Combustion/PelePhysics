#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


def test():
    import pyre.units

    parser = pyre.units.parser()
    parser = pyre.units.parser()
    
    print "The 7 base SI units:"
    print "             meter: %s" % parser.parse("meter")
    print "          kilogram: %s" % parser.parse("kilogram")
    print "            second: %s" % parser.parse("second")
    print "            ampere: %s" % parser.parse("ampere")
    print "            kelvin: %s" % parser.parse("kelvin")
    print "              mole: %s" % parser.parse("mole")
    print "           candela: %s" % parser.parse("candela")
    print
    print "The 22 SI derived units with special names:"
    print "            radian: %s" % parser.parse("radian")
    print "         steradian: %s" % parser.parse("steradian")
    print "             hertz: %s" % parser.parse("hertz")
                                                     
    print "            newton: %s" % parser.parse("newton")
    print "            pascal: %s" % parser.parse("pascal")
    print "             joule: %s" % parser.parse("joule")
    print "              watt: %s" % parser.parse("watt")
                                                     
    print "           coulomb: %s" % parser.parse("coulomb")
    print "              volt: %s" % parser.parse("volt")
    print "             farad: %s" % parser.parse("farad")
    print "               ohm: %s" % parser.parse("ohm")
    print "           siemens: %s" % parser.parse("siemens")
    print "             weber: %s" % parser.parse("weber")
    print "             tesla: %s" % parser.parse("tesla")
    print "             henry: %s" % parser.parse("henry")
                                                     
    print "    degree Celcius: %s" % parser.parse("celcius")
                                                     
    print "             lumen: %s" % parser.parse("lumen")
    print "               lux: %s" % parser.parse("lux")
                                                     
    print "         becquerel: %s" % parser.parse("becquerel")
    print "              gray: %s" % parser.parse("gray")
    print "           sievert: %s" % parser.parse("sievert")
    print "             katal: %s" % parser.parse("katal")

    return


# main
if __name__ == "__main__":
    test()


# version
__id__ = "$Id$"

# End of file 
