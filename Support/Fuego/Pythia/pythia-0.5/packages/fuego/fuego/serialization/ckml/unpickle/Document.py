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

from pyre.xml.DocumentNode import DocumentNode


class Document(DocumentNode):
    def mechanism(self):
        return self._mechanism

    def onKinetics(self, mechanism):
        self._mechanism = mechanism
        return

    def __init__(self, mechanism, source):
        DocumentNode.__init__(self, source, documentNodes())
        self._mechanism = mechanism

        return


# helpers


def documentNodes():

    from .Arrhenius import Arrhenius
    from .Atom import Atom
    from .Composition import Composition
    from .Duplicate import Duplicate
    from .Efficiencies import Efficiencies
    from .Element import Element
    from .Enhancement import Enhancement
    from .Kinetics import Kinetics
    from .LandauTeller import LandauTeller
    from .LowPressure import LowPressure
    from .NASA import NASA
    from .NASA_a1 import NASA_a1
    from .NASA_a2 import NASA_a2
    from .NASA_a3 import NASA_a3
    from .NASA_a4 import NASA_a4
    from .NASA_a5 import NASA_a5
    from .NASA_a6 import NASA_a6
    from .NASA_a7 import NASA_a7
    from .Product import Product
    from .Reactant import Reactant
    from .Reaction import Reaction
    from .ReactionRate import ReactionRate
    from .ReactionUnits import ReactionUnits
    from .Reagents import Reagents
    from .Reverse import Reverse
    from .ReverseLandauTeller import ReverseLandauTeller
    from .Species import Species
    from .SRI import SRI
    from .Thermo import Thermo
    from .TROE import TROE

    # primitives

    nodes = [
        Kinetics,
        Element,
        Species,
        Composition,
        Atom,
        Thermo,
        NASA,
        NASA_a1,
        NASA_a2,
        NASA_a3,
        NASA_a4,
        NASA_a5,
        NASA_a6,
        NASA_a7,
        Reaction,
        Duplicate,
        ReactionUnits,
        Reagents,
        Reactant,
        Product,
        ReactionRate,
        Efficiencies,
        Enhancement,
        Arrhenius,
        Reverse,
        LowPressure,
        SRI,
        TROE,
        LandauTeller,
        ReverseLandauTeller,
    ]

    return nodes


# version
__id__ = "$Id$"

# End of file
