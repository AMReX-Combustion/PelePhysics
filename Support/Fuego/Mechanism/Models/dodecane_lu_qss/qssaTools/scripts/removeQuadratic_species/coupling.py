def intersection(lst1, lst2):
    return list(set(lst1).intersection(lst2))

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def mergeSpeciesCoeff(speciesList,coeffList):
    # is there duplicate in species list
    speciesListUnique = list(set(speciesList))
    duplicate = (len(speciesList) != len(speciesListUnique))
    if duplicate:
       newCoeffList = []
       for species in speciesListUnique:
           indexDuplicates = list_duplicates_of(speciesList,species)
           newCoeffList.append(sum([coeffList[i] for i in indexDuplicates]))
       speciesList[:] = speciesListUnique
       coeffList[:] = newCoeffList

# Find possible problems
def identifyQSSCoupling(mechanism,qssSpeciesList):
   problematicQSS = []
   for ir, reaction in enumerate(mechanism.reaction()):
       reactant = [entry[0] for entry in reaction.reactants]
       coeffReac = [entry[1] for entry in reaction.reactants]
       mergeSpeciesCoeff(reactant,coeffReac)
       product = [entry[0] for entry in reaction.products]
       coeffProd = [entry[1] for entry in reaction.products]
       mergeSpeciesCoeff(product,coeffProd)
       # Check for quadratic relations
       speciesQSSAInvolved = intersection(reactant,qssSpeciesList)
       sumCoeffQSSAInvolved = sum([coeffReac[i] for i in [reactant.index(species) for species in speciesQSSAInvolved]])
       if sumCoeffQSSAInvolved>1:
           #print("ERROR LEFT " + reaction.equation())
           problematicQSS += speciesQSSAInvolved
       if reaction.reversible:
           speciesQSSAInvolved = intersection(product,qssSpeciesList)
           sumCoeffQSSAInvolved = sum([coeffProd[i] for i in [product.index(species) for species in speciesQSSAInvolved]])
           if sumCoeffQSSAInvolved>1:
               #print("ERROR RIGHT " + reaction.equation())
               problematicQSS += speciesQSSAInvolved

   problematicQSS = list(set(problematicQSS))
   return problematicQSS

# is there a problem
def QSSCoupling(mechanism,qssSpeciesList):
   error = 0
   for ir, reaction in enumerate(mechanism.reaction()):
       reactant = [entry[0] for entry in reaction.reactants]
       coeffReac = [entry[1] for entry in reaction.reactants]
       mergeSpeciesCoeff(reactant,coeffReac)
       product = [entry[0] for entry in reaction.products]
       coeffProd = [entry[1] for entry in reaction.products]
       mergeSpeciesCoeff(product,coeffProd)
       # Check for quadratic relations
       speciesQSSAInvolved = intersection(reactant,qssSpeciesList)
       sumCoeffQSSAInvolved = sum([coeffReac[i] for i in [reactant.index(species) for species in speciesQSSAInvolved]])
       if sumCoeffQSSAInvolved>1:
           error = 1
           return error
       if reaction.reversible:
           speciesQSSAInvolved = intersection(product,qssSpeciesList)
           sumCoeffQSSAInvolved = sum([coeffProd[i] for i in [product.index(species) for species in speciesQSSAInvolved]])
           if sumCoeffQSSAInvolved>1:
               error = 1
               return error

   return error

