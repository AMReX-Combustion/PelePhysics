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


from __future__ import print_function
def test():
    import pyre.units.SI as si
    
    print("The 7 base SI units:")
    print("             meter: %s" % si.meter)
    print("          kilogram: %s" % si.kilogram)
    print("            second: %s" % si.second)
    print("            ampere: %s" % si.ampere)
    print("            kelvin: %s" % si.kelvin)
    print("              mole: %s" % si.mole)
    print("           candela: %s" % si.candela)
    print()
    print("The 22 SI derived units with special names:")
    print("            radian: %s" % si.radian)
    print("         steradian: %s" % si.steradian)
    print("             hertz: %s" % si.hertz)
                                                     
    print("            newton: %s" % si.newton)
    print("            pascal: %s" % si.pascal)
    print("             joule: %s" % si.joule)
    print("              watt: %s" % si.watt)
                                                     
    print("           coulomb: %s" % si.coulomb)
    print("              volt: %s" % si.volt)
    print("             farad: %s" % si.farad)
    print("               ohm: %s" % si.ohm)
    print("           siemens: %s" % si.siemens)
    print("             weber: %s" % si.weber)
    print("             tesla: %s" % si.tesla)
    print("             henry: %s" % si.henry)
                                                     
    print("    degree Celcius: %s" % si.celcius)
                                                     
    print("             lumen: %s" % si.lumen)
    print("               lux: %s" % si.lux)
                                                     
    print("         becquerel: %s" % si.becquerel)
    print("              gray: %s" % si.gray)
    print("           sievert: %s" % si.sievert)
    print("             katal: %s" % si.katal)

    return


# main
if __name__ == "__main__":
    test()


# version
__id__ = "$Id$"

# End of file 
