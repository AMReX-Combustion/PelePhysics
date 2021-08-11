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

from builtins import object


class Species(object):
    def symbol(self):
        return self._symbol

    def molecularWeight(self):
        return self._molecularWeight

    def thermalCalculator(self):
        return self._thermalProperties

    def concentration(self):
        return self._concentration

    def updateConcentration(self, value):
        self._concentration = value

    def cp_R(self, T):
        return self._thermalProperties.cp_R(T.value)

    def h_RT(self, T):
        return self._thermalProperties.h_RT(T.value)

    def s_R(self, T):
        return self._thermalProperties.s_R(T.value)

    def g_RT(self, T):
        return self._thermalProperties.g_RT(T.value)

    def __init__(self, symbol, molecularWeight, thermalProperties):
        from pyre.handbook.units.SI import meter, mole

        self._symbol = symbol
        self._concentration = 0.0 * mole / meter ** 3
        self._molecularWeight = molecularWeight
        self._thermalProperties = thermalProperties

        return

    def __str__(self):
        return self._symbol


# version
__id__ = "$Id$"

#
# End of file
