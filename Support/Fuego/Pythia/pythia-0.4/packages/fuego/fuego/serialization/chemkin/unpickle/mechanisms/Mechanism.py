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


class Mechanism(object):


    from MechanismExceptions import DuplicateElement, DuplicateSpecies, DuplicateThermalProperties, DuplicateTransProperties


    # housekeeping
    
    def source(self):
        return self._source


    def printStatistics(self):
        print "Mechanism '%s'" % self._source
        print "    elements:", self._elements.size()
        print "     species:", self._species.size()
        print "      thermo:", self._thermoDb.size()
        print "      trans:", self._transDb.size()
        print "   reactions:", self._reactions.size()


    # elements

    def declareElement(self, element):
        symbol = element.symbol
        duplicate = self._elements.find(symbol)

        self._elements.element(element)

        if duplicate and element.locator:
            raise self.DuplicateElement(symbol)
        
        return


    def element(self, symbol=None):
        return self._elements.find(symbol)


    def elementDeclarations(self):
        return self._elements


    # species

    def declareSpecies(self, species):
        symbol = species.symbol
        duplicate = self._species.find(symbol)

        self._species.species(species)

        if duplicate:
            raise self.DuplicateSpecies(symbol)

        return


    def species(self, symbol=None):
        return self._species.find(symbol)


    def speciesDeclarations(self):
        return self._species


    # thermo


    def declareThermalProperties(self, species):
        symbol = species.symbol
        duplicate = self._thermoDb.find(symbol)

        self._thermoDb.species(species)

        if duplicate:
            raise self.DuplicateThermo(symbol)

        return


    def thermalProperties(self, species=None):
        prop = self._thermoDb.find(species)
        if not prop:
            return self._externalDb.find(species)


    def thermalDatabase(self):
        return self._thermoDb


    def thermoAll(self):
        return self._thermoDb.all(1)


    def thermoDone(self):
        if self._thermoDb.all():
            return

        self._externalDb = self._readExternalThermoDatabase()
        return


    def thermoRange(self, range=None):
        return self._thermoDb.range(range)


    # trans


    def declareTransProperties(self, species):
        symbol = species.symbol
        duplicate = self._transDb.find(symbol)

        self._transDb.species(species)

        if duplicate:
            raise self.DuplicateTrans(symbol)

        return


    def transProperties(self, species=None):
        prop = self._transDb.find(species)
        #if not prop:
        #    return self._externalDb.find(species)


    def transDatabase(self):
        return self._transDb


    def transAll(self):
        return self._transDb.all(1)


    def transDone(self):
        if self._transDb.all():
            return
        else:
            print "Someting went wrong with your transport"

        #self._externalDb = self._readExternalThermoDatabase()
        #return


    # reactions

    def declareReaction(self, reaction):
        self._reactions.reaction(reaction)
        return

    def reaction(self, species=None):
        return self._reactions.find(species)


    def reactionDeclarations(self):
        return self._reactions
        
        
    # other methods  

    def __init__(self, source):
        self._source = source

        self._externalDb = None

        from ElementDb import ElementDb
        self._elements = ElementDb()

        from SpeciesDb import SpeciesDb
        self._species = SpeciesDb()
        
        from ThermoDb import ThermoDb
        self._thermoDb = ThermoDb()

        from TransDb import TransDb
        self._transDb = TransDb()

        from ReactionDb import ReactionDb
        self._reactions = ReactionDb()

        from ElementDeclaration import ElementDeclaration
        self.elementFactory = ElementDeclaration
        
        from SpeciesDeclaration import SpeciesDeclaration
        self.speciesFactory = SpeciesDeclaration

        from ThermoDeclaration import ThermoDeclaration
        self.thermoFactory = ThermoDeclaration
        
        from ReactionDeclaration import ReactionDeclaration
        self.reactionFactory = ReactionDeclaration
        
        return


    def _readExternalThermoDatabase(self):
        import os
        import fuego
        from ExternalThermo import ExternalThermo

        therm = fuego.unpickle.chemkin.externalThermoDatabase()
        externalDb = open(therm, "r")

        mechanism = ExternalThermo(therm)

        pyre.diagnostics.category("chemkin-parser").maxHits(100000000)
        pyre.diagnostics.report(
            "chemkin-parser",
            "loading external thermo database '%s'" % externalDb.name)

        diagnosticState = pyre.diagnostics.category("chemkin-parser").active()

        #pyre.diagnostics.category("chemkin-parser").deactivate()
        mechanism = pyre.chemistry.unpickle.unpickleChemkin(externalDb, mechanism)
        if diagnosticState:
            pyre.diagnostics.category("chemkin-parser").activate()
            
        pyre.diagnostics.report(
            "chemkin-parser",
            "done loading external thermo database '%s'" % externalDb.name)

        db = mechanism.thermalDatabase()
        return db


# version
__id__ = "$Id$"

# End of file
