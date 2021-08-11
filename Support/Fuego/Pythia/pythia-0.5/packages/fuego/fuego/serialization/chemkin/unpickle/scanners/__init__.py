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


from __future__ import absolute_import


def sections():
    from .Sections import Sections

    return Sections()


def elements():
    from .Elements import Elements

    return Elements()


def species():
    from .Species import Species

    return Species()


def qss_species():
    from .QssSpecies import QssSpecies

    return QssSpecies()


def thermo():
    from .Thermo import Thermo

    return Thermo()


def trans():
    from .Trans import Trans

    return Trans()


def reactions():
    from .Reactions import Reactions

    return Reactions()


def parameters():
    from .Parameters import Parameters

    return Parameters()


# version
__id__ = "$Id$"

# End of file
