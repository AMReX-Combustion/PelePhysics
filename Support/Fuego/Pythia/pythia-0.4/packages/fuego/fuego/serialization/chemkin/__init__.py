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


def format():
    return "chemkin"


def extensions():
    return [".ck2"]


def parser():
    from unpickle.parsers.ChemkinParser import ChemkinParser
    return ChemkinParser()


def pickler():
    from pickle.ChemkinPickler import ChemkinPickler
    return ChemkinPickler()


# version
__id__ = "$Id$"

#  End of file 
