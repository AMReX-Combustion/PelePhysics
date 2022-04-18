import argparse
import os
import sys
sys.path.append(os.path.join(sys.path[0],'../utils'))
from QSSspecies import getListSpecies 
from coupling import *
from FMC import FMC
import fuego

# CLI
parser = argparse.ArgumentParser(description='Remove quadratic coupling according to method 2')
parser.add_argument('-mi', '--inputMechanism', type=str, metavar='', required=False, help='Mechanism to process', default='output/qssa_nostar.inp')
parser.add_argument('-th', '--thermoFile', type=str, metavar='', required=False, help='Thermodynamic file', default='output/therm_nostar.dat')
parser.add_argument('-tr', '--tranFile', type=str, metavar='', required=False, help='Transport file', default='output/tran_nostar.dat')
parser.add_argument('-nqss', '--nonQSSSpecies', type=str, metavar='', required=False, help='Filename with non QSS species names', default='output/non_qssa_list_nostar.txt')
parser.add_argument('-o', '--outputFolder', type=str, metavar='', required=False, help='Where to store by product of the preprocessing',default='output')
parser.add_argument('-mo', '--outputMechanism', type=str, metavar='', required=False, help='Mechanism to write', default='qssa_final.inp')
args = parser.parse_args()


# Set up filenames
input_mech_filename = args.inputMechanism
therm_filename = args.thermoFile
tran_filename = args.tranFile
non_qssa_species = args.nonQSSSpecies
outputFolder = args.outputFolder
output_mech_filename = args.outputMechanism

# Load skeletal mechanism
app = FMC()
app.inventory.mechanism = input_mech_filename
app.inventory.thermo = therm_filename
app.inventory.trans = tran_filename
mechanism = fuego.serialization.mechanism()
mechanism.externalThermoDatabase(app.inventory.thermo)
mechanism.externalTransDatabase(app.inventory.trans)
mechanism = fuego.serialization.load(filename=app.inventory.mechanism, format='chemkin', mechanism=mechanism)
reactionIndex = mechanism._sort_reactions()

# Get list of intended QSS species
species_non_qssa = getListSpecies(non_qssa_species)
species_all = getListSpecies(app.inventory.mechanism)
species_qssa = list(set(species_all) - set(species_non_qssa))


# Get reactions to cancel
toCancel = identifyReactionsToCancel(mechanism,species_qssa)
print(toCancel)

# Write new mechanism
f = open(input_mech_filename,'r')
lines = f.readlines()
f.close()

outputFolder = outputFolder
mech_final_filename = output_mech_filename
f = open(os.path.join(outputFolder,mech_final_filename),'w+')
line_iter = iter(lines)
for line in line_iter:
    line_write=line
    for ir in toCancel:
        if mechanism.reaction()[ir].equation().replace(' ','') in line:
            if 'f' in toCancel[ir] and 'r' in toCancel[ir]:
                #Dont print the reaction
                line_write = ''
                break
            elif 'f' in toCancel[ir]:
                if mechanism.reaction()[ir].reversible:
                    if not mechanism.reaction()[ir].low and not mechanism.reaction()[ir].thirdBody:
                        #Remove the entire reaction (this needs improvement)
                        line_write = ''
                        break
                    else:
                        line_write=''
                        #Dont write the next three lines
                        next(line_iter)
                        next(line_iter)
                        next(line_iter)
                        break
                else:
                    line_write=''
                    break
            elif 'r' in toCancel[ir]:
                    line_write = line.replace('<=>','=>')
                    print(mechanism.reaction()[ir].equation())
                    print(line_write)
                    break
            else:
                line_write=line
                break
    f.write(line_write)
f.close()




#Load new mech
app = FMC()
app.inventory.mechanism = os.path.join(outputFolder,mech_final_filename)
app.inventory.thermo = therm_filename
app.inventory.trans = tran_filename
mechanism = fuego.serialization.mechanism()
mechanism.externalThermoDatabase(app.inventory.thermo)
mechanism.externalTransDatabase(app.inventory.trans)
mechanism = fuego.serialization.load(filename=app.inventory.mechanism, format='chemkin', mechanism=mechanism)
reactionIndex = mechanism._sort_reactions()


if QSSCoupling(mechanism,species_qssa) == 0:
    print ('\n\n\nSUCCESS')
else:
    print ('\n\n\nFAILURE')

