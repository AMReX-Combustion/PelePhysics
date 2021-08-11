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

from weaver.mills.XMLMill import XMLMill


class CkmlPickler(XMLMill):
    def __init__(self):
        XMLMill.__init__(self)
        return

    def _renderDocument(self, mechanism, options=None):

        self._rep += ["", "<!DOCTYPE kinetics>", "", "<kinetics>", ""]

        self._indent()
        self._renderElements(mechanism)
        self._renderSpecies(mechanism)
        self._renderReactions(mechanism)
        self._outdent()

        self._rep += ["</kinetics>"]

        return

    def _renderElements(self, mechanism):
        self._write(self.line("elements"))
        for element in mechanism.element():
            symbol = element.symbol
            aw = element.weight

            tag = '<element id="%s"' % symbol
            if aw:
                tag += ' atomicWeight="%g"' % aw
            tag += "/>"
            self._write(tag)

        self._rep.append("")

        return

    def _renderSpecies(self, mechanism):
        self._write(self.line("species"))
        for species in mechanism.species():
            symbol = species.symbol
            phase = species.phase

            tag = '<species id="%s"' % symbol
            if phase:
                tag += ' phase="%s"' % phase
            tag += ">"
            self._write(tag)

            self._indent()

            composition = species.composition
            if composition:
                self._write("<composition>")

                self._indent()
                for element, coefficient in species.composition:
                    tag = '<atom element="%s"' % element
                    if coefficient != 1:
                        tag += ' coefficient="%d"' % coefficient
                    tag += "/>"
                    self._write(tag)

                self._outdent()
                self._write("</composition>")

            thermo = species.thermo
            if thermo:
                self._write("<thermo>")
                self._indent()
                for model in thermo:
                    self._write(
                        '<NASA lowT="%g" highT="%g">'
                        % (model.lowT, model.highT)
                    )
                    self._indent()

                    i = 1
                    for p in model.parameters:
                        self._write("<a%d>%g</a%d>" % (i, p, i))
                        i += 1

                    self._outdent()
                    self._write("</NASA>")

                self._outdent()
                self._write("</thermo>")

            self._outdent()
            self._write("</species>")

        self._rep.append("")

        return

    def _renderReactions(self, mechanism):
        self._write(self.line("reactions"))

        i = 0
        for reaction in mechanism.reaction():

            i += 1
            self._rep.append("")
            self._write(self.line(reaction.equation()))
            tag = '<reaction id="%d"' % i
            if reaction.reversible:
                tag += ' reversible="true"'
            if reaction.thirdBody:
                species, coefficient = reaction.thirdBody
                if species == "<mixture>":
                    species = "mixture"
                tag += ' thirdBody="%s"' % species
            if reaction.falloff:
                tag += ' falloff="true"'

            tag += ">"
            self._write(tag)

            self._indent()

            if reaction.duplicate:
                self._write("<duplicate/>")

            if reaction.units:
                tag = "<units "

                activation = reaction.units.get("activation")
                if activation:
                    tag += ' activation="%s"' % activation

                prefactor = reaction.units.get("prefactor")
                if prefactor:
                    tag += ' prefactor="%s"' % prefactor

                tag += "/>"

                self._write(tag)

            self._write("<reagents>")
            self._indent()
            for species, coefficient in reaction.reactants:
                tag = '<reactant species="%s"' % species
                if coefficient != 1:
                    tag += ' coefficient="%d"' % coefficient
                tag += "/>"
                self._write(tag)

            for species, coefficient in reaction.products:
                tag = '<product species="%s"' % species
                if coefficient != 1:
                    tag += ' coefficient="%d"' % coefficient
                tag += "/>"
                self._write(tag)

            self._outdent()
            self._write("</reagents>")

            self._write("<rate>")
            self._indent()

            arrhenius = reaction.arrhenius
            if arrhenius:
                tag = '<arrhenius A="%g" beta="%g" E="%g"/>' % tuple(arrhenius)
                self._write(tag)

            low = reaction.low
            if low:
                tag = '<low A="%g" beta="%g" E="%g"/>' % tuple(low)
                self._write(tag)

            reverse = reaction.rev
            if reverse:
                tag = '<reverse A="%g" beta="%g" E="%g"/>' % tuple(reverse)
                self._write(tag)

            sri = reaction.sri
            if sri:
                param = tuple(sri)
                if len(param) == 3:
                    tag = '<sri a="%g" b="%g" c="%g"/>' % tuple(sri)
                elif len(param) == 5:
                    tag = '<sri a="%g" b="%g" c="%g" d="%g" e="%g"/>' % tuple(
                        sri
                    )
                else:
                    import pyre

                    pyre.debug.Firewall.hit(
                        "poorly formed SRI parameter tuple"
                    )
                self._write(tag)

            troe = reaction.troe
            if troe:
                tag = '<troe a="%g" T3s="%g" Ts="%g"' % tuple(troe[:3])
                if len(troe) == 4:
                    tag += ' T2s="%g"' % troe[3]
                tag += "/>"
                self._write(tag)

            lt = reaction.lt
            if lt:
                tag = '<lt B="%g" C="%g"/>' % tuple(lt)
                self._write(tag)

            rlt = reaction.rlt
            if rlt:
                tag = '<rlt B="%g" C="%g"/>' % tuple(rlt)
                self._write(tag)

            efficiencies = reaction.efficiencies
            if efficiencies:
                self._write("<efficiencies>")
                self._indent()
                for species, coefficient in reaction.efficiencies:
                    tag = '<enhancement species="%s" factor="%g"/>' % (
                        species,
                        coefficient,
                    )
                    self._write(tag)

                self._outdent()
                self._write("</efficiencies>")

            self._outdent()
            self._write("</rate>")

            self._outdent()
            self._write("</reaction>")

        self._rep.append("")
        return


# version
__id__ = "$Id$"

#  End of file
