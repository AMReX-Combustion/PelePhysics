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


from __future__ import absolute_import
def format():
    return "native"


def extensions():
    return [".py"]


def parser():
    from .NativeParser import NativeParser
    return NativeParser()


def pickler():
    from .NativePickler import NativePickler
    return NativePickler()


# version
__id__ = "$Id$"

#  End of file 
