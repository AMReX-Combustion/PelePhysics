import numpy as np
import argparse
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


# CLI
parser = argparse.ArgumentParser(description='Remove star from species names')
parser.add_argument('-sk', '--skeletalMechanism', type=str, metavar='', required=True, help='Skeletal mechanism')
parser.add_argument('-m', '--mechanismToProcess', type=str, metavar='', required=False, help='Mechanism file to process', default='output/qssa.inp')
parser.add_argument('-th', '--thermoToProcess', type=str, metavar='', required=True, help='Thermo file to process')
parser.add_argument('-tr', '--tranToProcess', type=str, metavar='', required=True, help='Transport file to process')
parser.add_argument('-nqss', '--nonQSSSpecies', type=str, metavar='', required=True, help='Filename with non QSS species names')
parser.add_argument('-sknos', '--skeletalMechanismNoStar', type=str, metavar='', required=False, help='Skeletal mechanism', default='skeletal_nostar.inp')
parser.add_argument('-mnos', '--mechanismToProcessNoStar', type=str, metavar='', required=False, help='Mechanism file to write', default='qssa_nostar.inp')
parser.add_argument('-thnos', '--thermoToProcessNoStar', type=str, metavar='', required=False, help='Thermo file to write', default='therm_nostar.dat')
parser.add_argument('-trnos', '--tranToProcessNoStar', type=str, metavar='', required=False, help='Transport file to write', default='tran_nostar.dat')
parser.add_argument('-nqssnos', '--nonQSSSpeciesNoStar', type=str, metavar='', required=False, help='Filename with non QSS species names and no star', default='non_qssa_list_nostar.txt')
parser.add_argument('-o', '--outputFolder', type=str, metavar='', required=False, help='Where to store by product of the preprocessing', default='output')
args = parser.parse_args()

# Set up filenames
skeletal_mech_filename = args.skeletalMechanism
mech_filename = args.mechanismToProcess
therm_filename = args.thermoToProcess
tran_filename = args.tranToProcess
non_qssa_list_filename = args.nonQSSSpecies
outputFolder = args.outputFolder
skeletal_mech_nostar_filename = os.path.join(outputFolder,args.skeletalMechanismNoStar)
mech_nostar_filename = os.path.join(outputFolder,args.mechanismToProcessNoStar)
therm_nostar_filename = os.path.join(outputFolder,args.thermoToProcessNoStar)
tran_nostar_filename = os.path.join(outputFolder,args.tranToProcessNoStar)
non_qssa_list_nostar_filename = os.path.join(outputFolder,args.nonQSSSpeciesNoStar)

# Get species
species_all = getListSpecies(mech_filename)

# Which species contain a star
starSpecies = []
for species in species_all:
    if "*" in species:
        starSpecies.append(species)


# Replace star in species name
# Skeletal Mech
write_nostar_inp(inputFile=skeletal_mech_filename,
                 starSpecies=starSpecies,
                 outputFile=skeletal_mech_nostar_filename)
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

# Non qssa
write_nostar_inp(inputFile=non_qssa_list_filename,
                 starSpecies=starSpecies,
                 outputFile=non_qssa_list_nostar_filename)
