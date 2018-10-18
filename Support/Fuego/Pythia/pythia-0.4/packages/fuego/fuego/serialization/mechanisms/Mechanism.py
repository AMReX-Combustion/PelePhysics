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

import journal


class Mechanism(object):


    from MechanismExceptions import DuplicateElement, DuplicateSpecies, DuplicateThermalProperties


    # housekeeping
    
    def name(self):
        return self._name


    def externalThermoDatabase(self, filename):
        self._thermdb = filename
        return

    def externalTransDatabase(self, filename):
        self._transdb = filename
        return


    def printStatistics(self):
        print "Mechanism '%s'" % self._name
        print "    elements:", self._elements.size()
        print "     species:", self._species.size()
        print "   reactions:", self._reactions.size()


    # elements

    def newElement(self, symbol, weight=None, locator=None):
        duplicate = self._elements.find(symbol)

        element = self._elements.element(symbol, weight, locator)

        if duplicate and element.locator:
            raise self.DuplicateElement(symbol)
        
        return element


    def element(self, symbol=None):
        return self._elements.find(symbol)


    # species

    def newSpecies(self, symbol, locator=None):
        duplicate = self._species.find(symbol)

        species = self._species.species(symbol, locator)

        if duplicate:
            raise self.DuplicateSpecies(symbol)

        return species


    def species(self, symbol=None):
        return self._species.find(symbol)


    # thermal properties are recorded directly in the species


    def thermoAll(self, flag=None):
        if not flag:
            self._externalDb = self._readExternalThermoDb()
            self._thermoRange = self._externalDb.thermoRange()
            
        return self._externalDb


    # trigger the ingestion of therm.dat
    def thermoDone(self):

        unresolvedSpecies = []

        for species in self._species.find():
            if not species.thermo:
                if not self._externalDb:
                    self._externalDb = self._readExternalThermoDb()

                resolution = self._externalDb.species(species.symbol)
                resolution.trans = species.trans

                if not resolution:
                    unresolvedSpecies.append(species)
                else:
                    self._info.log(
                        "resolving species '%s' against '%s'" % (species.symbol, self._thermdb))
                    self._species.replace(species.symbol, species, resolution)

        if unresolvedSpecies:
            warning = journal.warning("fuego")
            warning.line("unresolved species in mechanism")
            warning.line("species: %s" % [ x.symbol for x in unresolvedSpecies])
                
        return 0


    def thermoRange(self, range=None):
        if range:
            self._thermoRange = range
        return self._thermoRange


    # reactions

    def newReaction(self, id, locator=None):
        return self._reactions.reaction(id, locator)


    def reaction(self, species=None, id=None):
        if not self._sorted:
            print '*** WARNING: reactions have not been sorted'
        return self._reactions.find(species, id)


    def _sort_reactions(self):
        n = [0]
        rs = []
        rs_unsorted = self._reactions.find()
        i = 0
        # troe
        for r in rs_unsorted:
            if r not in rs:
                if r.low and r.troe and not r.rev: 
                    i+=1
                    r.orig_id = r.id
                    r.id = i
                    rs.append(r)
        n.append(i)
        # sri
        for r in rs_unsorted:
            if r not in rs:
                if r.low and r.sri and not r.rev: 
                    i+=1
                    r.orig_id = r.id
                    r.id = i
                    rs.append(r)
        n.append(i)
        # lindemann
        for r in rs_unsorted:
            if r not in rs:
                if r.low and not r.rev: 
                    i+=1
                    r.orig_id = r.id
                    r.id = i
                    rs.append(r)
        n.append(i)
        # three-body:
        for r in rs_unsorted:
            if r not in rs:
                if r.thirdBody and not r.low and not r.rev: 
                    i+=1
                    r.orig_id = r.id
                    r.id = i
                    rs.append(r)
        n.append(i)
        # simplest case
        for r in rs_unsorted:
            if r not in rs:
                if not r.rev and not r.low and not r.thirdBody: 
                    i+=1
                    r.orig_id = r.id
                    r.id = i
                    rs.append(r)
        n.append(i)
        # everything else
        for r in rs_unsorted:
            if r not in rs:
                i+=1
                r.orig_id = r.id
                r.id = i
                rs.append(r)
        n.append(i)

        for r in rs:
            self._reactions.replace2(r,r.id-1,r)
        self._sorted = True

        return n


    # other methods  

    def __init__(self, name=""):
        from ElementSet import ElementSet
        from SpeciesSet import SpeciesSet
        from ReactionSet import ReactionSet

        self._name = name
        self._thermdb = "therm.dat"
        self._transdb = "tran.dat"
        self._externalDb = None

        self._elements = ElementSet()
        self._species = SpeciesSet()
        self._reactions = ReactionSet()

        self._thermoRange = ()

        self._info = journal.debug("fuego.serialization")

        self._sorted = False
        
        return


    # swallow an external thermo database
    def _readExternalThermoDb(self):
        import fuego
        filename = self._thermdb
        db = fuego.serialization.loadThermoDatabase(filename, format="chemkin")
        return db


    def dump(self):
        print
        print "Statistics:"
        print "-----------"
        self.printStatistics()

        print
        print "Elements:"
        print "---------"
        for element in self.element():
            print "%6s: %s" % (element.symbol, element)
            
        print
        print "Species:"
        print "---------"
        for species in self.species():
            print "%10s: %s" % (species.symbol, species)

        print
        print "Reactions:"
        print "----------"
        i = 1
        for reaction in self.reaction():
            print "%4d: %s" % (i, reaction)
            i += 1

        return




# version
__id__ = "$Id$"

# End of file
