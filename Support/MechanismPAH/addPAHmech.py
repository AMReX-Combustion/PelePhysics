

"""
This script compares an input Chemkin mechanism, transport, and thermodynamic data
with a PAH sub-mechanism and creates a YAML file of the combined mechanism
Care is taken to ensure duplicate species and reactions are not created
@author: Landon Owen
"""

import logging
import cleanfiles as clf
from getdata import *
import compareenthalpy
import compareIDT as cidt
import argparse
import numpy as np
import re
import os
import sys
import matplotlib.pyplot as plt
import cantera.ck2yaml as cky
import cantera as ct

def getspecnames(mechfile):
    speciesnames = []
    with open(mechfile) as infile:
        for line in infile:
            if (len(line) > 0):
                if (line[0:7] == "SPECIES"):
                    spline = next(infile)
                    while (not ("END" in spline)):
                        specs = spline.split()
                        for sp in specs:
                            speciesnames.append(sp)
                        spline = next(infile)
                    break
    return speciesnames
def checkyaml(mechfile, newthermo, newtransport,permissive):
    cky_final = cky.Parser()
    if(permissive):
        cky_final.warning_as_error = False
    cky_final.load_chemkin_file(mechfile)
    cky_final.load_chemkin_file(newthermo)
    cky_final.load_chemkin_file(newtransport)
    numMechreacts = len(cky_final.reactions)
    yamlfile = os.path.splitext(mechfile)[0] + '.yaml'
    surf = cky_final.write_yaml('gas', yamlfile)
    logger = logging.getLogger(__name__)
    dupreacts = []
    dupeq = []
    try:
        logger.info('Validating mechanism...')
        gas = ct.Solution(yamlfile)
        logger.info('PASSED')
        return [False, dupreacts, dupeq]
    except RuntimeError as e:
        logger.info('FAILED')
        msg = str(e)
        if 'Undeclared duplicate reactions' in msg:
            for line in msg.split('\n'):
                match = re.match('>.*# Reaction ([0-9]+)', line)
                if match:
                    cindx = int(match.group(1))-1
                    pahindx = cindx - numMechreacts
                    dupreacts.append(pahindx)
                    eq = line.split(":")[-1].split("#")[0]
                    dupeq.append(eq)
                    print("Undeclated duplicate reaction: " + str(eq))
        else:
            logger.warning(e)
    return [True, dupreacts, dupeq]

parser = argparse.ArgumentParser()

parser.add_argument("--plotH", help="Flag to plot enthalpy values", type=bool, default=False)
parser.add_argument("--permissive", help="Flag to allow some recoverable parsing errors", type=bool, default=False)
parser.add_argument("--mech", help="Mechanism input file", type=str, required=True)
parser.add_argument("--thermo", help="Thermodynamic input file", type=str, required=True)
parser.add_argument("--transport", help="Transport input file", type=str, required=True)
parser.add_argument("--fuelname", help="Name of primary fuel used in mechanism", type=str, required=True)
parser.add_argument("--newmech", help="Mechanism output file name", type=str, default="mechanism.inp")
parser.add_argument("--newthermo", help="Thermodynamic output file name", type=str, default="thermo.dat")
parser.add_argument("--newtransport", help="Transport output file name", type=str, default="trans.dat")

args = parser.parse_args()
mechspecs = getspecnames(args.mech)
# First create clean files that don't have extra species
clf.cleanfiles(mechspecs, args.thermo, args.transport, args.newthermo, args.newtransport)
pahpath = "PAHmechfiles"
pahthermo = pahpath + '/PAH.chthermo' # This should remain the same
pahmech = pahpath + '/PAH.chmech' # This should remain the same
pahtrans = pahpath + '/PAH.chtrans' # This should remain the same
pahspecs = getspecnames(pahmech)
# Compare the enthalpy values of various species to see if they are in fact the same
# This returns a list of species names that must be changed from the PAH chemistry
[modnamePAH, modnameMech] = compareenthalpy.compareenthalpy(pahthermo, args.newthermo, args.plotH)
# This is a list of species that are entirely new to the existing mechanism and must be added
newPAHspec = []
for pahspec in pahspecs:
    if (pahspec in modnamePAH):
        indx = modnamePAH.index(pahspec)
        print("Must rename " + pahspec + " in PAH module to " + modnameMech[indx])
    elif (not (pahspec in mechspecs)):
        newPAHspec.append(pahspec)
        print("Must add species " + str(pahspec))
# Add thermo data for new species
addnewthermdata(pahthermo, args.newthermo, newPAHspec)

# Now add transport data for new species
addnewtransdata(pahtrans, args.newtransport, newPAHspec)
[modpahreacts, fullmechfile] = updatemechfile(args.mech, args.newmech, pahmech, newPAHspec, modnamePAH, modnameMech)
# Use Cantera to import the newly created mechanism and check for duplicate reactions
removed = 0
[found_error, dupreact, dupeq] = checkyaml(args.newmech, args.newthermo, args.newtransport, args.permissive)
while (found_error):
    removed += 1
    pahindx = max(dupreact)
    print("Removing equation " + dupeq[dupreact.index(pahindx)]  + " from PAH reactions")
    modpahreacts.pop(pahindx)
    with open(args.newmech, 'w') as infile:
        for line in fullmechfile:
            infile.write(line)
        for line in modpahreacts:
            infile.write(line)
        infile.write("END\n")
    [found_error, dupreact, numMechreacts] = checkyaml(args.newmech, args.newthermo, args.newtransport, args.permissive)
print("Removed " + str(removed) + " reactions due to duplication")


newyamlfile = os.path.splitext(args.newmech)[0] + '.yaml'
# Create a YAML from original mechanism
[found_error, dupreact, numMechreacts] = checkyaml(args.mech, args.thermo, args.transport, args.permissive)
origyamlfile = os.path.splitext(args.mech)[0] + '.yaml'
# Test the ignition delay time at 1 atm and 20 atm
# Output values into a plot file

numatm = 1
cidt.compareIDT(origyamlfile, newyamlfile, numatm, args.fuelname)
numatm = 20
cidt.compareIDT(origyamlfile, newyamlfile, numatm, args.fuelname)
