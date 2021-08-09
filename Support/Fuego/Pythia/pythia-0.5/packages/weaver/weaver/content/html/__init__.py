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
def weaver():
    from .DocumentWeaver import DocumentWeaver
    return DocumentWeaver()


def body():
    from .StandardBody import StandardBody
    return StandardBody()


def section():
    from .Section import Section
    return Section()


def stateSelector():
    from .StatesSelection import StatesSelection
    return StatesSelection()


def textBlock(message):
    from .TextBlock import TextBlock
    return TextBlock(message)


# version
__id__ = "$Id$"

#  End of file 
