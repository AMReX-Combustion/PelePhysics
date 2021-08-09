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

from weaver.mills.PythonMill import PythonMill


class NativePickler(PythonMill):


    def __init__(self):
        PythonMill.__init__(self)
        return


    def _renderDocument(self, mechanism, options=None):
        self._rep += [
            '',
            'from fuego.serialization.mechanisms.NASA import NASA',
            '',
            ''
            ]

        self._renderElements(mechanism)
        self._renderSpecies(mechanism)
        self._renderReactions(mechanism)
        return


    def _end(self):
        self._timestamp()
        self._rep += self.footer()
        return


    def _renderElements(self, mechanism):
        self._write()
        self._write(self.line(self.separator()))
        self._write(self.line('elements'))
        for element in mechanism.element():
            symbol = element.symbol
            aw = element.weight

            line = 'mechanism.newElement("%s"' % symbol
            if aw:
                line += ', %s' % aw
            line += ')'

            self._write(line)

        return


    def _renderSpecies(self, mechanism):

        self._write('')
        self._write(self.line(self.separator()))
        self._write(self.line('species'))

        for species in mechanism.species():
            symbol = species.symbol
            self._write()
            self._write(self.line('%s' % symbol))

            line = '_species = mechanism.newSpecies("%s")' % symbol
            self._write(line)

            line = '_species.phase = "%s"' % species.phase
            self._write(line)

            line = '_species.composition = %s' % species.composition
            self._write(line)

            for thermo in species.thermo:
                line = "_thermo = NASA(%g, %g)" % (thermo.lowT, thermo.highT)
                self._write(line)
                line = "_thermo.parameters = %r" % (thermo.parameters,)
                self._write(line)
                line = "_species.thermo.append(_thermo)"
                self._write(line)
                
            self._outdent()

        return


    def _renderReactions(self, mechanism):

        self._write('')
        self._write(self.line(self.separator()))
        self._write(self.line('reactions'))

        for reaction in mechanism.reaction():
            self._write()
            self._write(self.line("%03d: %s" % (reaction.id, reaction.equation())))

            line = '_reaction = mechanism.newReaction(%d)' % reaction.id
            self._write(line)

            line = '_reaction.reactants = %s' % reaction.reactants
            self._write(line)
            line = '_reaction.products = %s' % reaction.products
            self._write(line)
            line = '_reaction.units = %s' % reaction.units
            self._write(line)

            if reaction.duplicate:
                line = '_reaction.duplicate = %s' % reaction.duplicate
                self._write(line)
            
            if reaction.reversible:
                line = '_reaction.reversible = %s' % reaction.reversible
                self._write(line)
            
            if reaction.thirdBody:
                line = '_reaction.thirdBody = %s' % repr(reaction.thirdBody)
                self._write(line)
            
            if reaction.falloff:
                line = '_reaction.falloff = %s' % reaction.falloff
                self._write(line)
            
            if reaction.efficiencies:
                line = '_reaction.efficiencies = %s' % reaction.efficiencies
                self._write(line)

            line = '_reaction.arrhenius = %s' % repr(reaction.arrhenius)
            self._write(line)
            
            if reaction.low:
                line = '_reaction.low = %s' % repr(reaction.low)
                self._write(line)

            if reaction.sri:
                line = '_reaction.sri = %s' % repr(reaction.sri)
                self._write(line)

            if reaction.troe:
                line = '_reaction.troe = %s' % repr(reaction.troe)
                self._write(line)

            if reaction.lt:
                line = '_reaction.lt = %s' % repr(reaction.lt)
                self._write(line)

            if reaction.rev:
                line = '_reaction.rev = %s' % repr(reaction.rev)
                self._write(line)

            if reaction.rlt:
                line = '_reaction.rlt = %s' % repr(reaction.rlt)
                self._write(line)

        return


# version
__id__ = "$Id$"

#  End of file 
