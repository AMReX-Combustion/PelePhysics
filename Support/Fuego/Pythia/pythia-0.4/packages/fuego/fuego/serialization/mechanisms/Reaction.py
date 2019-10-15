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

from Entity import Entity


class Reaction(Entity):


    def equation(self):
        str = self._reagents(self.reactants)

        if self.reversible:
            str += " <=> "
        else:
            str += " => "

        str += self._reagents(self.products)
        return str


    def __init__(self, id, locator=None):
        Entity.__init__(self, id, locator)

        self.products = []
        self.reactants = []
        self.efficiencies = []
        self.ford = []

        self.arrhenius = None
        self.units = {}

        self.lt = None
        self.rev = None
        self.rlt = None

        self.duplicate = None
        self.thirdBody = None
        self.reversible = None

        self.falloff = None
        self.low = None
        self.sri = None
        self.troe = None

        self.radiation = None

        return


    def _reagents(self, composition):
        terms = []
        for species, factor in composition:
            str = ""
            if factor != 1:
                str += "%f " % factor
            str += species

            terms.append(str)

        reaction = " + ".join(terms)

        if self.thirdBody:
            thirdBody = ""
            species, factor = self.thirdBody
            if species == '<mixture>':
                species = 'M'

            if self.falloff:
                thirdBody += ' (+'
            else:
                thirdBody += ' + '

            if factor != 1:
                thirdBody += "%f" % factor

            thirdBody += species

            if self.falloff:
                thirdBody += ')'
        
            reaction += thirdBody

        return reaction



    def __str__(self):

        str = self.equation()
        str += ", id=%d" % self.id
        str += ", source=" + Entity.__str__(self)
        return str


# version
__id__ = "$Id$"

# End of file
