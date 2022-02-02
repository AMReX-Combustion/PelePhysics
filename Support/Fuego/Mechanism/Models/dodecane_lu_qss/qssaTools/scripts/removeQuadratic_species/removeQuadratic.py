from FMC import FMC
import fuego
from QSSspecies import getListSpecies 
from coupling import *
import sys
sys.path.append('scripts/utils')
import myparser
import os
import numpy as np

def writeQSSA_inp(skeletal_inp_filename,species_qssa,outputFolder,qssa_inp_filename):
    f = open(skeletal_inp_filename,'r+')
    lines_skeletal = f.readlines()
    f.close()

    os.makedirs(outputFolder,exist_ok=True)
    f = open(outputFolder+'/'+qssa_inp_filename,'w+')
    # Write elements
    startWrite = False
    for iline, line in enumerate(lines_skeletal):
        if line.startswith('ELEMENT'):
            startWrite = True
        if startWrite:
            f.write(line)
        if startWrite and line.startswith('END'):
            break
    f.write('\n')
    # Write species
    startWrite = False
    for iline, line in enumerate(lines_skeletal):
        if line.startswith('SPECIES'):
            startWrite = True
        if startWrite:
            f.write(line)
        if startWrite and line.startswith('END'):
            break
    f.write('\n')
    # Write QSS
    f.write('QSS')
    f.write('\n')
    string = ''
    for ispec, species in enumerate(species_qssa):
        #print(ispec)
        if not (ispec+1)%4==0 and (ispec+1)<len(species_qssa):
           string += species
           for i in range(17-len(species)):
               string += ' '
        else:
           string += species
           for i in range(17-len(species)):
               string += ' '
           string += '\n'
           f.write(string)
           #print(string)
           string=''
    f.write('END')
    f.write('\n')
    # Write the rest
    for line in lines_skeletal[iline+1:]:
        f.write(line)
    f.close()

# Parse input
inpt = myparser.parseInputFile()

# Set up filenames
skeletal_mech_filename = inpt['skeletal_mech']
mech_filename = inpt['mech']
therm_filename = inpt['therm']
tran_filename = inpt['tran']


# Load skeletal mechanism
app = FMC()
app.inventory.mechanism = mech_filename
app.inventory.thermo = therm_filename
app.inventory.trans = tran_filename
mechanism = fuego.serialization.mechanism()
mechanism.externalThermoDatabase(app.inventory.thermo)
mechanism.externalTransDatabase(app.inventory.trans)
mechanism = fuego.serialization.load(filename=app.inventory.mechanism, format='chemkin', mechanism=mechanism)
reactionIndex = mechanism._sort_reactions()

# Get list of intended QSS species
species_non_qssa = getListSpecies(inpt['non_qssa_list'])
species_all = getListSpecies(app.inventory.mechanism)
species_qssa = list(set(species_all) - set(species_non_qssa))

# Get possible issues
problematic_species_qssa = identifyQSSCoupling(mechanism,species_qssa)

# Find all possible combinaisons of species to remove from QSS list to get rid of quadratic relations
print('Here are possible qssa species to remove : ')
import itertools
tryspecies_qssa=[]  
qssa_remove_proposal = []
for L in range(0, len(problematic_species_qssa)+1):
    for subset in itertools.combinations(problematic_species_qssa, L):
        # Make a new list of qssa species
        tryspecies_qssa[:] = species_qssa     
        for speciesRemove in list(subset):
            tryspecies_qssa.remove(speciesRemove)
        # Check if still creates problem
        if QSSCoupling(mechanism,tryspecies_qssa)==0:
            # This combinaison works
            # Does this combinaison contain entirely another successful combinaison
            # If yes, then we are removing too many species, do not include it as solution
            add = True
            for success_qssa_found in qssa_remove_proposal:
                if set(success_qssa_found).issubset(list(subset)):
                    add=False
                    break
            if add:
                qssa_remove_proposal.append(list(subset))
                print(subset)


# Where the smallest set to remove
len_species_remove = [len(subset) for subset in qssa_remove_proposal]
ind=np.argwhere(np.array(len_species_remove)==np.amin(np.array(len_species_remove)))[0][0]
qssa_species_remove = qssa_remove_proposal[ind] 

# Make new set of qss species
newSpecies_qssa = species_qssa.copy()
for species in qssa_species_remove:
    newSpecies_qssa.remove(species)


# Write new mechanism
outputFolder = inpt['outputFolder']
mech_final_filename = inpt['mech_final']
writeQSSA_inp(skeletal_mech_filename,newSpecies_qssa,outputFolder,mech_final_filename)

