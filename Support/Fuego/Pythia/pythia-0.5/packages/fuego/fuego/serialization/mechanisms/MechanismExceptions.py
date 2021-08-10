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

# elements

class DuplicateElement(Exception):

    def __init__(self, symbol):
        self._symbol = symbol
        return


    def __str__(self):
        return "duplicate element '%s'" % self._symbol


# species

class DuplicateSpecies(Exception):


    def __init__(self, symbol):
        self._symbol = symbol
        return


    def __str__(self):
        return "duplicate species '%s'" % self._symbol


# qss species

class DuplicateQssSpecies(Exception):


    def __init__(self, symbol):
        self._symbol = symbol
        return


    def __str__(self):
        return "duplicate QSS species '%s'" % self._symbol


# thermo

class DuplicateThermalProperties(Exception):


    def __init__(self, symbol):
        self._symbol = symbol
        return


    def __str__(self):
        return "duplicate thermodynamic properties for '%s'" % self._symbol


# version
__id__ = "$Id$"

#
# End of file
