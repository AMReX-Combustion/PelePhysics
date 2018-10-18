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


def sections():
    from Sections import Sections
    return Sections()


def elements():
    from Elements import Elements
    return Elements()


def species():
    from Species import Species
    return Species()


def thermo():
    from Thermo import Thermo
    return Thermo()


def trans():
    from Trans import Trans
    return Trans()


def reactions():
    from Reactions import Reactions
    return Reactions()


def parameters():
    from Parameters import Parameters
    return Parameters()


# version
__id__ = "$Id$"

# End of file
