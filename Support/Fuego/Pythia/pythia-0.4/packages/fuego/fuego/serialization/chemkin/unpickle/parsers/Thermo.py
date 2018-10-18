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

from BaseParser import BaseParser


class Thermo(BaseParser):


    # the interesting tokens

    def aTemperatureRange(self, token):
        if not self._thermoAll:
            msg = "Unexpected temperature range definitions without THERMO ALL"
            self.onWarning(msg, self.locator())

        self._range = token.range
        self._mechanism.thermoRange(self._range)

        return 0


    def aThermoLine(self, token):
        if not self._range and self._thermoAll:
            if not self._thermoAllWarned:
                msg = "THERMO ALL: Expected temperature range definition, not species info"
                self.onWarning(msg, self.locator())
                self._thermoAllWarned = 1
            
        # is this the next valid thermo line in the sequence?
        id = token.id - 1
        if id != self._nextId:
            msg = "Unexpected thermo line: found %d while expecting %d" % \
                  (token.id, self._nextId + 1)
            self.onError(msg, self.locator())

        # dispatch to appropriate line info parser
        self._lineParsers[id](token) 
        
        # next valid thermo line is...
        self._nextId = (self._nextId + 1) % 4

        return 0


    # transitions

    def aThermoSection(self, token):
        self._info.log("thermo parser: section start")

        self._thermoAll = token._all
        self._mechanism.thermoAll(self._thermoAll)
        self._range = self._mechanism.thermoRange()
        self._parse(self._scanner, self._tokenizer)

        return 0


    def anEndSection(self, token):
        BaseParser.anEndSection(self, token)

        if self._thermoAll:
            species = self._mechanism.species()

            candidates = []
            for s in species:
                if not s.thermo:
                    candidates.append(s.symbol)

            if candidates:
                msg = "no species information for the folowing species: %s" % candidates
                self.onWarning(msg, self.locator())

        return 1


    def onEndOfFile(self):
        # if not self._thermoAll:
            # msg = "this mechanism requires an external thermo database"
            # self.onWarning(msg)

        self._mechanism.thermoDone()
            
        return 1


    # other methods

    def __init__(self, mechanism, tokenizer):
        import pyre
        import fuego
        BaseParser.__init__(self, mechanism)

        self._tokenizer = tokenizer
        self._scanner = fuego.serialization.chemkin.unpickle.scanners.thermo()

        # Private data
        self._range = ()
        self._thermoAll = 0
        self._parameters = []
        self._currentRange = None

        import re
        from fuego.serialization.chemkin.unpickle.tokens.RegularExpressions import species
        self._speciesScanner = re.compile(species)

        self._nextId = 0 #line ids are zero-based

        self._currentSpecies = None

        self._lineParsers = [
            self._parseLine1, self._parseLine2, self._parseLine3, self._parseLine4
            ]

        self._thermoAllWarned = 0

        return


    def _parseLine1(self, token):
        text = token.text

        match = self._speciesScanner.match(text[0:18])
        if not match:
            msg = "Could not match a valid species name in '%s'" % text[0:18]
            self.onError(msg, self.locator())

        speciesName = match.group()

        species = self._mechanism.species(speciesName)
        if not species:
            msg = "thermo section: undeclared species '%s'" % speciesName
            self.onWarning(msg, self.locator())
            species = self._mechanism.newSpecies(speciesName)
            
        species.locator(self.locator())

        # Parse the element coefficient in columns 24-43 (zero-based)
        for i in range(0, 4):
            offset = 24+i*5
            self._extractComposition(token, text, offset, species.composition)

        # Get the phase
        phase = text[44].lower()
        if phase not in ["s", "l", "g"]:
            msg = "Unkown phase code '%s'" % phase
            locator = self.locator()
            locator.column = 44
            self.onError(msg, self.locator())

        species.phase = phase

        # Get the temperature intervals
        lowT = self._extractFloat(token, text, 45, 10)
        highT = self._extractFloat(token, text, 55, 10)
        midT = self._extractFloat(token, text, 65, 10, optional=1)
        if midT == None:
            midT = self._range[1]

        self._currentRange = (lowT, midT, highT)

        # The extra possible element,coef pair
        self._extractComposition(token, text, 73, species.composition)

        # Save this information
        self._currentSpecies = species
        
        return


    def _parseLine2(self, token):
        if not self._currentSpecies: return

        text = token.text

        # extract the high T range parametrization
        self._parameters = []
        for i in range(0, 5):
            number = self._extractFloat(token, text, i*15, 15)
            self._parameters.append(number)
        
        return


    def _parseLine3(self, token):
        if not self._currentSpecies: return

        text = token.text

        # finish extracting the high T range parametrization
        for i in range(0, 2):
            number = self._extractFloat(token, text, i*15, 15)
            self._parameters.append(number)

        # store in the species
        self._currentSpecies.thermalParametrization(
            "NASA", self._currentRange[1], self._currentRange[2], self.locator(),
            self._parameters
            )
        
        # extract the first part of the low T parameters
        self._parameters = []
        for i in range(2, 5):
            number = self._extractFloat(token, text, i*15, 15)
            self._parameters.append(number)
        
        return


    def _parseLine4(self, token):
        species = self._currentSpecies
        if not species: return

        text = token.text

        # finish extracting the low T range parametrization
        for i in range(0, 4):
            number = self._extractFloat(token, text, i*15, 15)
            self._parameters.append(number)
        
        # store in the species
        self._currentSpecies.thermalParametrization(
            "NASA", self._currentRange[0], self._currentRange[1], self.locator(),
            self._parameters
            )
        
        return
            

    def _extractFloat(self, token, text, offset, width, optional=0):
        str = text[offset:offset+width].strip()
        if not str:
            if optional: return None
            msg = "Expected a required numeric field instead of '%s'" % text[offset:offset+width]
            locator = self.locator()
            locator.column = offset
            self.onError(msg, locator)
            

        try:
            value = float(str)
        except ValueError:
            msg = "Could not convert '%s' into a number" % text[offset:offset+width]
            locator = self.locator()
            locator.column = offset
            self.onError(msg, locator)

        return value
            

    def _extractComposition(self, token, text, offset, composition):

        # Extract the coefficient first:
        #    some files have junk in the element slot when the coefficient is 0
        coef = text[offset+2:offset+5].strip()

        # The coefficient could be blank, which means 0
        if not coef: return

        coef = self._extractFloat(token, text, offset+2, 3)

        # Extract the element name
        name = text[offset:offset+2].strip()
        if name and coef:
            element = self._mechanism.element(name)
            if not element:
                msg = "Element '%s' not declared in an ELEMENT section" % name
                locator = self.locator()
                locator.column = offset
                self.onWarning(msg, locator)

            composition.append( (name, coef) )

        return
        

# version
__id__ = "$Id$"

#
# End of file
