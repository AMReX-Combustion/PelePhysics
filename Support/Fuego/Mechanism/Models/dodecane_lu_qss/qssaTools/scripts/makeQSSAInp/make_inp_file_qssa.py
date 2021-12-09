import numpy as np
import sys
sys.path.append('scripts/utils')
import myparser
import os

def getListSpecies(filename):
    f = open(filename,'r+')
    lines = f.readlines()
    f.close()
    # Get species names    
    species = []
    startRead = False
    for line in lines:
        if line.startswith('SPECIES'):
           startRead = True
        if line.startswith('END') and startRead:
           break
        if not line.startswith('SPECIES') and startRead: 
           species += line.split()
    return species

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
skeletal_inp_filename = inpt['skeletal_inp']
species_non_qssa_filename = inpt['species_non_qssa']
outputFolder = inpt['outputFolder'] 
qssa_inp_filename = inpt['qssa_inp']

# Get species
species_non_qssa = getListSpecies(species_non_qssa_filename)
species_all = getListSpecies(skeletal_inp_filename)
species_qssa = list(set(species_all) - set(species_non_qssa))
#print(species_qssa)

# Write QSSA mech
writeQSSA_inp(skeletal_inp_filename=skeletal_inp_filename,
              species_qssa=species_qssa,
              outputFolder=outputFolder,
              qssa_inp_filename=qssa_inp_filename)

