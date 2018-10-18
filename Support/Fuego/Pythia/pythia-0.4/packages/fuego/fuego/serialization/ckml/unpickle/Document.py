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

    from Kinetics import Kinetics

    from Element import Element

    from Species import Species
    from Composition import Composition
    from Atom import Atom
    from Thermo import Thermo
    from NASA import NASA
    from NASA_a1 import NASA_a1
    from NASA_a2 import NASA_a2
    from NASA_a3 import NASA_a3
    from NASA_a4 import NASA_a4
    from NASA_a5 import NASA_a5
    from NASA_a6 import NASA_a6
    from NASA_a7 import NASA_a7

    from Reaction import Reaction
    from Duplicate import Duplicate
    from ReactionUnits import ReactionUnits
    from Reagents import Reagents
    from Reactant import Reactant
    from Product import Product
    from ReactionRate import ReactionRate
    from Efficiencies import Efficiencies
    from Enhancement import Enhancement
    from Arrhenius import Arrhenius
    from Reverse import Reverse
    from LowPressure import LowPressure
    from SRI import SRI
    from TROE import TROE
    from LandauTeller import LandauTeller
    from ReverseLandauTeller import ReverseLandauTeller

    # primitives
    
    nodes = [
        Kinetics,

        Element,

        Species, Composition, Atom, Thermo,
        NASA, NASA_a1, NASA_a2, NASA_a3, NASA_a4, NASA_a5, NASA_a6, NASA_a7,

        Reaction, Duplicate, ReactionUnits, Reagents, Reactant, Product,
        ReactionRate, Efficiencies, Enhancement,
        Arrhenius, Reverse, LowPressure, SRI, TROE, LandauTeller, ReverseLandauTeller,
        ]

    return nodes


# version
__id__ = "$Id$"

# End of file
