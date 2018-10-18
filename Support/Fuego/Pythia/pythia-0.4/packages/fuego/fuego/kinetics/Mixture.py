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

import pyre
from EntitySet import EntitySet


class Mixture(EntitySet):


    def moles(self, symbol, value):
        
        species = self.find(symbol)
        species.updateConcentration(value)
        
        mixture = self.find("<mixture>")
        mixture.updateConcentration(mixture.concentration() + value)

        return


    def __init__(self, mechanism):
        EntitySet.__init__(self)

        from ElementSet import ElementSet
        from Species import Species

        elements = ElementSet(mechanism.element())

        for proxy in mechanism.species():

            symbol = proxy.symbol
            weight = molecularWeight(elements, proxy.composition)
            thermo = thermalProperties(proxy.thermo)
            
            species = Species(symbol, weight, thermo)
            self.insert(proxy.symbol, species)

        mixture = Species("<mixture>", 0, None)
        self.insert(mixture.symbol(), mixture)

        return


# helpers

def molecularWeight(elements, composition):
    weight = 0.0
        
    for symbol, coefficient in composition:
        weight += coefficient * elements.find(symbol).atomicWeight()
       
    return weight


def thermalProperties(thermo):
    # for now assume that thermo is a list of two NASA parametrizations

    r1, r2 = tuple(thermo)

    if r1.lowT < r2.lowT:
        low = r1
        high = r2
    else:
        low = r2
        high = r1

    pyre.debug.Firewall.verify(low.highT == high.lowT, "bad thermal properties")
    range = (low.lowT, low.highT, high.highT)
    nasa = pyre.models.thermodynamics.nasa(range, low.parameters, high.parameters)
    
    return nasa


# version
__id__ = "$Id$"

# End of file
