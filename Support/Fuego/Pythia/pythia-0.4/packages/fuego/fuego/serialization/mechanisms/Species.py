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

from Entity import Entity


class Species(Entity):


    def thermalParametrization(self, type, lowT, highT, locator, parameters):
        if type == "NASA":
            from NASA import NASA
            model = NASA(lowT, highT, locator)
            model.parameters = parameters
        else:
            import pyre
            pyre.debug.Firewall.hit("unknown thermal parametrization type '%s'" % type)
            return

        self.thermo.append(model)
        return

    def transParametrization(self, type, EPS, SIG, DIP, POL, ZROT, locator, parameters):
        #if type == "0":
        from TRANLIN import TRANLIN
        model = TRANLIN(EPS, SIG, DIP, POL, ZROT, locator)
        model.parameters = parameters
        #elif type == "1":
        #    model = TRANNL(EPS, SIG, DIP, POL, ZROT, locator)
        #    model.parameters = parameters
        #elif type == "2":
        #    model = TRANNL(EPS, SIG, DIP, POL, ZROT, locator)
        #    model.parameters = parameters
        #else:
        #    import pyre
        #    pyre.debug.Firewall.hit("unknown transport parametrization type '%s'" % type)
        #    return

        self.trans.append(model)
        return

    def __init__(self, id, symbol, locator=None):
        Entity.__init__(self, id, locator)
        self.symbol = symbol
        self.phase = None
        self.composition = []
        self.thermo = []
        self.trans = []
        return


    def __str__(self):
        str = "symbol='%s'" % self.symbol
        str += ", phase='%s'" % self.phase
        str += ", composition=%s" % self.composition


        for p in self.thermo:
            str += ", thermo=([%g, %g]: %s)" % (p.lowT, p.highT, p.parameters)
        for p in self.trans:
            str += ", trans=([%s, %s, %s, %s, %s]: %s)" % ((p.eps), p.sig, p.dip, p.pol, p.zrot, p.parameters)

        str += ", source=" + Entity.__str__(self)
        return str


# version
__id__ = "$Id$"

# End of file
