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

    #HARI======================================
    #   These are routines to sort reactions
    #   and improve memory locality
    #==========================================
    def _reorder_reaction_set_tsp(self,rset):

        import numpy as npy

        nReactions = len(rset)
        nSpecies = len(self.species())

        reactionmat=npy.zeros((nReactions,nSpecies))

        for i,reaction in zip(range(nReactions),rset):

            agents = list(set(reaction.reactants+reaction.products))

            for a in agents:
                symbol, coefficient = a
                reactionmat[i][self.species(symbol).id]=coefficient
        
        new_to_old_map=self._tsp_solve(reactionmat,0.001)
        #new_to_old_map=self._cluster_solve(reactionmat)

        print(new_to_old_map)

        return(new_to_old_map)

    def _sort_reactions_within_type_tsp(self,n):

        #note--------------------------------------------------------
        #Species ids, ie sp.id starts with 0
        #while reaction ids, ie reaction.id starts with 1
        #although when doing mechanism.reaction(id=i), i should be
        #index starting with 0. this is because the "find()" function 
        #in reactionSet class queries within the entity array
        #------------------------------------------------------------

        #sort within each type
        #=====================
        rs    = self._reactions.find()

        itroe      = n[0:2]
        isri       = n[1:3]
        ilindemann = n[2:4]
        i3body     = n[3:5] 
        isimple    = n[4:6]
        ispecial   = n[5:7]

        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        troe_order    = self._reorder_reaction_set_tsp(rs[itroe[0]:itroe[1]])+itroe[0]
        sri_order     = self._reorder_reaction_set_tsp(rs[isri[0]:isri[1]])+isri[0]
        lind_order    = self._reorder_reaction_set_tsp(rs[ilindemann[0]:ilindemann[1]])+ilindemann[0]
        thbody_order  = self._reorder_reaction_set_tsp(rs[i3body[0]:i3body[1]])+i3body[0]
        simple_order  = self._reorder_reaction_set_tsp(rs[isimple[0]:isimple[1]])+isimple[0]
        special_order = self._reorder_reaction_set_tsp(rs[ispecial[0]:ispecial[1]])+ispecial[0]

        new_to_old_map = troe_order.tolist()+sri_order.tolist()+\
        lind_order.tolist()+thbody_order.tolist()+simple_order.tolist()+special_order.tolist()

        self._reorder_reactions_from_map(new_to_old_map)

    def _sort_reactions_within_type_random(self,n):

        import numpy as npy
        #note--------------------------------------------------------
        #Species ids, ie sp.id starts with 0
        #while reaction ids, ie reaction.id starts with 1
        #although when doing mechanism.reaction(id=i), i should be
        #index starting with 0. this is because the "find()" function 
        #in reactionSet class queries within the entity array
        #------------------------------------------------------------

        #sort within each type
        #=====================
        rs    = self._reactions.find()

        itroe      = n[0:2]
        isri       = n[1:3]
        ilindemann = n[2:4]
        i3body     = n[3:5] 
        isimple    = n[4:6]
        ispecial   = n[5:7]

        ntroe      = itroe[1]      - itroe[0]
        nsri       = isri[1]       - isri[0]
        nlindemann = ilindemann[1] - ilindemann[0]
        n3body     = i3body[1]     - i3body[0]
        nsimple    = isimple[1]    - isimple[0]
        nspecial   = ispecial[1]   - ispecial[0]

        troe_order    = npy.random.permutation(ntroe)      + itroe[0]
        sri_order     = npy.random.permutation(nsri)       + isri[0]
        lind_order    = npy.random.permutation(nlindemann) + ilindemann[0]
        thbody_order  = npy.random.permutation(n3body)     + i3body[0]
        simple_order  = npy.random.permutation(nsimple)    + isimple[0]
        special_order = npy.random.permutation(nspecial)   + ispecial[0]

        new_to_old_map = troe_order.tolist()+sri_order.tolist()+\
        lind_order.tolist()+thbody_order.tolist()+simple_order.tolist()+special_order.tolist()

        self._reorder_reactions_from_map(new_to_old_map)

    def _reorder_reactions_from_map(self,new_to_old_map):

        rs    = self._reactions.find()
        rsnew = []

        for i in range(len(new_to_old_map)):
            r = rs[new_to_old_map[i]]
            r.id = i+1 #id should start with 1
            rsnew.append(r)
        
        for r in rsnew:
            self._reactions.replace2(r,r.id-1,r)

    def _reorder_species_from_map(self,new_to_old_map):

        from SpeciesSet import SpeciesSet
        import copy

        nSpecies = len(self.species())
        spnew=SpeciesSet()

        #reorder species
        for i in range(nSpecies):
            for sp in self.species():
                if(sp.id == new_to_old_map[i]):
                    break

            sp_temp=copy.deepcopy(sp)
            sp_temp.id=i
            spnew.insert(sp_temp.symbol, sp_temp)

        self._species=spnew


    def _get_reaction_matrix(self):

        import numpy as npy

        nSpecies = len(self.species())
        nReactions = len(self.reaction())

        reactionmat=npy.zeros((nReactions,nSpecies))

        for i in range(nReactions):

            reaction = self.reaction(id=i) #here id has to start from 0
            agents = list(set(reaction.reactants+reaction.products))
            efficiencies = reaction.efficiencies

            for a in agents:
                symbol, coefficient = a
                reactionmat[i][self.species(symbol).id]=coefficient
            
            for ii, eff in enumerate(efficiencies):
                symbol, efficiency = eff
                reactionmat[i][self.species(symbol).id]=1.0

        return(reactionmat)

    def _cluster_solve(self,mat):

        from sklearn.cluster import AgglomerativeClustering
        import numpy as npy

        new_to_old_map=npy.array([])

        if(mat.shape[0] > 1):

            nclus=mat.shape[0]/4
            #nclus=2

            clustering = AgglomerativeClustering(n_clusters=nclus, compute_full_tree=True, affinity='l1', linkage='average')
            y=clustering.fit_predict(mat)

            for i in range(nclus):
                for j in range(len(y)):
                    if(y[j]==i):
                        new_to_old_map = npy.append(new_to_old_map,j)

            new_to_old_map=new_to_old_map.astype(int)


        else:
            new_to_old_map=npy.arange(mat.shape[0])

        return(new_to_old_map)


    def _tsp_solve(self,mat,improvement_threshold):

        import numpy as npy

        #===============================================================
        # Calculate the euclidian distance in n-space of the route r traversing cities c, ending at the path start.
        path_distance = lambda r,c: npy.sum([npy.linalg.norm(c[r[p]]-c[r[p-1]],1) for p in range(len(r))])
        # Reverse the order of all elements from element i to element k in array r.
        two_opt_swap = lambda r,i,k: npy.concatenate((r[0:i],r[k:-len(r)+i-1:-1],r[k+1:len(r)]))

        def two_opt(cities,improvement_threshold): # 2-opt Algorithm adapted from https://en.wikipedia.org/wiki/2-opt

            route = npy.arange(cities.shape[0]) # Make an array of row numbers corresponding to cities.
            improvement_factor = 1 # Initialize the improvement factor.
            best_distance = path_distance(route,cities) # Calculate the distance of the initial path.
            while improvement_factor > improvement_threshold: # If the route is still improving, keep going!
                distance_to_beat = best_distance # Record the distance at the beginning of the loop.
                for swap_first in range(1,len(route)-2): # From each city except the first and last,
                    for swap_last in range(swap_first+1,len(route)): # to each of the cities following,
                        new_route = two_opt_swap(route,swap_first,swap_last) # try reversing the order of these cities
                        new_distance = path_distance(new_route,cities) # and check the total distance with this modification.
                        if new_distance < best_distance: # If the path distance is an improvement,
                            route = new_route # make this the accepted best route
                            best_distance = new_distance # and update the distance corresponding to this route.
                improvement_factor = 1 - best_distance/distance_to_beat # Calculate how much the route has improved.
            return route # When the route is no longer improving substantially, stop searching and return the route.
        #===============================================================

        if(len(mat) > 0):
            nrows=mat.shape[0]
            ncols=mat.shape[1]
            newmat=npy.zeros((nrows+1,ncols))
            newmat[1:(nrows+1),:]=mat
            order=two_opt(newmat,improvement_threshold)
            return(order[1:(nrows+1)]-1)
        else:
            return(npy.array([]))


    def _sort_species_ids_tsp(self):

        import numpy as npy

        rmat=self._get_reaction_matrix()
        new_to_old_map=self._tsp_solve(npy.transpose(rmat),0.001)
        #new_to_old_map=self._cluster_solve(npy.transpose(rmat))

        self._reorder_species_from_map(new_to_old_map)

    def _sort_species_ids_random(self):

        import numpy as npy

        nSpecies = len(self.species())
        rmat=self._get_reaction_matrix()
        new_to_old_map = npy.random.permutation(nSpecies)

        self._reorder_species_from_map(new_to_old_map)

    #===================================================================


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
