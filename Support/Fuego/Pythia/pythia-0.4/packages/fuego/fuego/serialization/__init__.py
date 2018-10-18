#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

from pickle import save, pickler, picklers
from unpickle import load, loadThermoDatabase, unpickler, unpicklers


def formats():
    from pickle import registrar as picklers
    from unpickle import registrar as unpicklers

    # return (input, output) registered formats
    return (unpicklers().registered(), picklers().registered())


def mechanism(filename=""):
    from mechanisms.Mechanism import Mechanism
    return Mechanism(filename)


# version
__id__ = "$Id$"

#  End of file 
