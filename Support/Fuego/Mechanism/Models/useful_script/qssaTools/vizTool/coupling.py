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


# Find possible problems
def identifyReactionsToCancel(mechanism,qssSpeciesList):
   problematicQSS = []
   toCancel = {}
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
          if not ir in toCancel:
              toCancel[ir] = []
          toCancel[ir].append('f')
       if reaction.reversible:
           speciesQSSAInvolved = intersection(product,qssSpeciesList)
           sumCoeffQSSAInvolved = sum([coeffProd[i] for i in [product.index(species) for species in speciesQSSAInvolved]])
           if sumCoeffQSSAInvolved>1:
               if not ir in toCancel:
                  toCancel[ir] = []
               toCancel[ir].append('r')

   return toCancel 


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


# Build Quad graph
def getQuadGraphLinks(mechanism,qssSpeciesList):
   links = []
   quadraticReac = [] 
   problematicQSSSpecies = []
   for ir, reaction in enumerate(mechanism.reaction()):
       reactant = [entry[0] for entry in reaction.reactants]
       coeffReac = [entry[1] for entry in reaction.reactants]
       mergeSpeciesCoeff(reactant,coeffReac)
       product = [entry[0] for entry in reaction.products]
       coeffProd = [entry[1] for entry in reaction.products]
       mergeSpeciesCoeff(product,coeffProd)
       # Identify QSS species involved
       reacQSSAInvolved = intersection(reactant,qssSpeciesList)
       sumCoeffReacQSSAInvolved = sum([coeffReac[i] for i in [reactant.index(species) for species in reacQSSAInvolved]])
       prodQSSAInvolved = intersection(product,qssSpeciesList)
       sumCoeffProdQSSAInvolved = sum([coeffProd[i] for i in [product.index(species) for species in prodQSSAInvolved]])
       if (prodQSSAInvolved) and (reacQSSAInvolved):
           for prod in prodQSSAInvolved:
               links.append((prod,'F'+str(ir)))
               for reac in reacQSSAInvolved:
                   links.append(('F'+str(ir),reac))
           if reaction.reversible:
               for reac in reacQSSAInvolved:
                   links.append((reac,'R'+str(ir)))
                   for prod in prodQSSAInvolved:
                       links.append(('R'+str(ir),prod))
       if (not prodQSSAInvolved) and (reacQSSAInvolved) and sumCoeffReacQSSAInvolved>1:
           for reac in reacQSSAInvolved:
               links.append(('F'+str(ir),reac))
       if (prodQSSAInvolved) and (not reacQSSAInvolved) and sumCoeffProdQSSAInvolved>1 and reaction.reversible:
           for prod in prodQSSAInvolved:
               links.append(('R'+str(ir),prod))

       if sumCoeffReacQSSAInvolved>1:
           # This is a quadratic coupling
           quadraticReac.append('F'+str(ir))
           problematicQSSSpecies += reacQSSAInvolved
           problematicQSSSpecies += prodQSSAInvolved

       if reaction.reversible:
           if sumCoeffProdQSSAInvolved>1:
               quadraticReac.append('R'+str(ir))
               problematicQSSSpecies += reacQSSAInvolved
               problematicQSSSpecies += prodQSSAInvolved
        

   return links, quadraticReac, problematicQSSSpecies


# Build directed graph
def getDiGraphLinks(mechanism,qssSpeciesList):
   links = []
   quadraticReac = [] 
   problematicQSSSpecies = []
   for ir, reaction in enumerate(mechanism.reaction()):
       reactant = [entry[0] for entry in reaction.reactants]
       coeffReac = [entry[1] for entry in reaction.reactants]
       mergeSpeciesCoeff(reactant,coeffReac)
       product = [entry[0] for entry in reaction.products]
       coeffProd = [entry[1] for entry in reaction.products]
       mergeSpeciesCoeff(product,coeffProd)
       # Identify QSS species involved
       reacQSSAInvolved = intersection(reactant,qssSpeciesList)
       prodQSSAInvolved = intersection(product,qssSpeciesList)
       sumCoeffReacQSSAInvolved = sum([coeffReac[i] for i in [reactant.index(species) for species in reacQSSAInvolved]])
       sumCoeffProdQSSAInvolved = sum([coeffProd[i] for i in [product.index(species) for species in prodQSSAInvolved]])
       if (prodQSSAInvolved) and (reacQSSAInvolved):
           for prod in prodQSSAInvolved:
               for reac in reacQSSAInvolved:
                   links.append((prod,reac))
           if reaction.reversible:
               for reac in reacQSSAInvolved:
                   for prod in prodQSSAInvolved:
                       links.append((reac,prod))
       if len(reacQSSAInvolved)>1:
           for reac1 in reacQSSAInvolved:
               for reac2 in reacQSSAInvolved:
                   if not reac1==reac2:
                       links.append((reac1,reac2))

       if len(prodQSSAInvolved)>1:
           for prod1 in prodQSSAInvolved:
               for prod2 in prodQSSAInvolved:
                   if not prod1==prod2:
                       links.append((prod1,prod2))


   return links

