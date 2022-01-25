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


def replace_all(origString,starSpecies):
    # Make a dict of how origin and new string
    replaceRuleDict = {}
    for species in starSpecies:
        replaceRuleDict[species] = species.replace("*","D")
    
    for i, j in replaceRuleDict.items():
        origString = origString.replace(i, j)
    return origString

def write_nostar_inp(inputFile,starSpecies,outputFile):
    
    with open(inputFile, 'r') as read_file:
        content = read_file.readlines()

    newName = species.replace("*","D")
    with open(outputFile, 'w') as write_file:
        for line in content:
            newLine = replace_all(line,starSpecies)
            write_file.write(newLine)
  
    return

# Parse input

inpt = myparser.parseInputFile()

# Set up filenames
mech_filename = inpt['mech']
therm_filename = inpt['therm']
tran_filename = inpt['tran']
outputFolder = inpt['outputFolder'] 
mech_nostar_filename = os.path.join(outputFolder,inpt['mech_nostar'])
therm_nostar_filename = os.path.join(outputFolder,inpt['therm_nostar'])
tran_nostar_filename = os.path.join(outputFolder,inpt['tran_nostar'])

# Get species
species_all = getListSpecies(mech_filename)

# Which species contain a star
starSpecies = []
for species in species_all:
    if "*" in species:
        starSpecies.append(species)


# Replace star in species name
# Mech
write_nostar_inp(inputFile=mech_filename,
                 starSpecies=starSpecies,
                 outputFile=mech_nostar_filename)
# Therm
write_nostar_inp(inputFile=therm_filename,
                 starSpecies=starSpecies,
                 outputFile=therm_nostar_filename)
# Tran
write_nostar_inp(inputFile=tran_filename,
                 starSpecies=starSpecies,
                 outputFile=tran_nostar_filename)
