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

class Element:


    def symbol(self):
        return self._symbol


    def atomicWeight(self):
        return self._atomicWeight


    def __init__(self, symbol, atomicWeight=None):
        # normalize symbol according to the conventions for PeriodicTable
        self._symbol = symbol.capitalize()

        if not atomicWeight:
            import pyre
            element = pyre.handbook.periodicTable().symbol(self._symbol)
            if not element:
                import journal
                journal.error("fuego").line("element '%s': no atomic weight information" % symbol)
                journal.error("fuego").log(
                    "attempted to look it up in the periodic table as '%s'" % self._symbol)
                raise pyre.exceptions.UnrecoverableError("unknown element '%s'" % symbol)
                
            else:
                atomicWeight = element.atomicWeight

        self._atomicWeight = atomicWeight
        
        return


    def __str__(self):
        return self._symbol


# version
__id__ = "$Id$"

#
# End of file
