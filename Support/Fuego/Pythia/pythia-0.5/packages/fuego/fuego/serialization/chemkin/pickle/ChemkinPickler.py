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

from weaver.mills.LineMill import LineMill


class ChemkinPickler(LineMill):


    names = ["chemkin"]


    def _renderDocument(self, mechanism, options=None):

        self.pickleElementSection(mechanism)
        self.pickleSpeciesSection(mechanism)
        self.pickleThermoSection(mechanism)
        self.pickleReactionSection(mechanism)

        return


    def pickleElementSection(self, mechanism):
        self._rep += [
            '',
            '! Element section',
            '',
            'Elements'
            ]

        line = " "*4
        
        for element in mechanism.element():
            symbol = element.symbol
            if len(line) + len(symbol) > 75:
                self._rep.append(line)
                line = " "*4

            line += " " + symbol

        self._rep.append(line)
        self._rep.append('End')

        return


    def pickleSpeciesSection(self, mechanism):
        self._rep += [
            '',
            '! Species section',
            '',
            'Species'
            ]

        line = " "*4
        for species in mechanism.species():
            symbol = species.symbol
            if len(line) + len(symbol) > 75:
                self._rep.append(line)
                line = " "*4

            line += " " + symbol

        self._rep.append(line)
        self._rep.append('End')

        return


    def pickleThermoSection(self, mechanism):
        self._rep += [
            '',
            '! Thermo section',
            '']

        line = 'Thermo'
        if mechanism.thermoAll():
            line += " All"
        self._rep.append(line)

        if mechanism.thermoRange():
            line = "%15.8g "*3 % mechanism.thermoRange()
            self._rep.append(line)

        format = "%15.8e"*5 + "%5d"

        for species in mechanism.species():

            if not species.thermo:
                continue

            self._rep.append("!")

            # compute line 1

            line_1 = "%-18s" % species.symbol + " "*6 

            composition = [
                "%-2s%3d" % (element, factor)
                for element, factor in species.composition]

            line_1 += "".join(composition[:min(len(composition), 4)])
            line_1 += (" "*5)*(max(0, 4-len(composition)))
            line_1 += species.phase.upper()

            line_1 += "%10.3f" % species.thermo[1].lowT
            line_1 += "%10.3f" % species.thermo[0].highT

            if species.thermo[1].highT != species.thermo[0].lowT:
                import journal
                journal.firewall("fuego").hit("bad mechanism")
                continue
                
            if species.thermo[1].lowT:
                line_1 += "%10.3f" % species.thermo[1].lowT
            else:
                line_1 += " "*10

            if len(composition) >= 5:
                line_1 += "%-2s%2d" % composition[4]
            else:
                line_1 += " "*4

            line_1 += "1"
            self._rep.append(line_1)

            # get the thermo parametrization

            highParameters = species.thermo[0].parameters
            lowParameters = species.thermo[1].parameters

            # compute line 2

            line_2 = ""
            line_2 += "%15.8e" % highParameters[0]
            line_2 += "%15.8e" % highParameters[1]
            line_2 += "%15.8e" % highParameters[2]
            line_2 += "%15.8e" % highParameters[3]
            line_2 += "%15.8e" % highParameters[4]
            line_2 += " "*4 + "2"

            self._rep.append(line_2)

            # compute line 3

            line_3 = ""
            line_3 += "%15.8e" % highParameters[5]
            line_3 += "%15.8e" % highParameters[6]
            line_3 += "%15.8e" % lowParameters[0]
            line_3 += "%15.8e" % lowParameters[1]
            line_3 += "%15.8e" % lowParameters[2]
            line_3 += " "*4 + "3"

            self._rep.append(line_3)

            # compute line 4
            line_4 = ""
            line_4 += "%15.8e" % lowParameters[3]
            line_4 += "%15.8e" % lowParameters[4]
            line_4 += "%15.8e" % lowParameters[5]
            line_4 += "%15.8e" % lowParameters[6]
            line_4 += " "*15
            line_4 += " "*4 + "4"

            self._rep.append(line_4)

        self._rep.append('')
        self._rep.append('End')

        return


    def pickleReactionSection(self, mechanism):
        self._rep.append('')
        self._rep.append('! Reaction section')
        self._rep.append('')
        self._rep.append('Reactions')

        i = 0

        for reaction in mechanism.reaction():
            i += 1
            self.pickleReaction(reaction, i)


        self._rep.append('')
        self._rep.append('End')

        return


    def pickleReaction(self, reaction, i):
        lines = []
        form = _printReagents(reaction, reaction.reactants)

        if reaction.reversible:
            form += " <=> "
        else:
            form += " => "

        form += _printReagents(reaction, reaction.products)

        line = "%-40s" % form

        line += "%10.3g" % reaction.arrhenius[0]
        line += "%10.3g" % reaction.arrhenius[1]
        line += "%10.3g" % reaction.arrhenius[2]
        line += " "*5 + "! %5d" % i
        lines.append(line)

        if reaction.efficiencies:
            efficiencies = "    "
            for species, coefficient in reaction.efficiencies:
                efficiencies += "%s / %4.2f / " % (species, coefficient + 1) # remember adjustment

            lines.append(efficiencies)

        if reaction.low:
            low = "    LOW /%s/" % _printParameters(reaction.low)
            lines.append(low)

        if reaction.troe:
            troe = "    TROE /%s/" % _printParameters(reaction.troe)
            lines.append(troe)

        if reaction.sri:
            sri = "    SRI /%s/" % _printParameters(reaction.sri)
            lines.append(sri)

        if reaction.rev:
            rev = "    REV /%s/" % _printParameters(reaction.rev)
            lines.append(rev)

        if reaction.lt:
            lt = "    LT /%s/" % _printParameters(reaction.lt)
            lines.append(lt)

        if reaction.rlt:
            rlt = "    RLT /%s/" % _printParameters(reaction.rlt)
            lines.append(rlt)

        if reaction.radiation:
            radiation = "    HV / %g /" % reaction.radiation
            lines.append(radiation)

        if reaction.duplicate:
            duplicate = "    DUPLICATE"
            lines.append(duplicate)

        self._rep += lines

        return lines


    def __init__(self, options=None):
        LineMill.__init__(self, '!', _FIRSTLINE)
        return


# helpers

_FIRSTLINE = '! -*- chemkin -*-'


def _printReagents(reaction, composition):
    terms = []
    for species, factor in composition:
        str = ""
        if factor != 1:
            str += "%d " % factor

        str += species

        terms.append(str)
        
    line = " + ".join(terms)

    if reaction.thirdBody:
        species, factor = reaction.thirdBody
        if species == '<mixture>':
            species = 'M'

        if reaction.falloff:
            line += ' (+'
        else:
            line += ' + '

        if factor != 1:
            line += "%d" % factor

        line += species

        if reaction.falloff:
            line += ')'
        
    return line


def _printParameters(ptuple):
    format = "%10.3e " * len(ptuple)
    return format % ptuple


# version
__id__ = "$Id$"

# End of file
