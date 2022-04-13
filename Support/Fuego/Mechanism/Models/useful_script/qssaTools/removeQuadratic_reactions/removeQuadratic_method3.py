from FMC import FMC
import fuego
import argparse
from QSSspecies import getListSpecies 
from coupling import *
import os

# CLI
parser = argparse.ArgumentParser(description='Remove quadratic coupling according to method 3')
parser.add_argument('-mi', '--inputMechanism', type=str, metavar='', required=False, help='Mechanism to process', default='output/qssa_nostar.inp')
parser.add_argument('-th', '--thermoFile', type=str, metavar='', required=False, help='Thermodynamic file', default='output/therm_nostar.dat')
parser.add_argument('-tr', '--tranFile', type=str, metavar='', required=False, help='Transport file', default='output/tran_nostar.dat')
parser.add_argument('-nqss', '--nonQSSSpecies', type=str, metavar='', required=False, help='Filename with non QSS species names', default='output/non_qssa_list_nostar.txt')
parser.add_argument('-o', '--outputFolder', type=str, metavar='', required=False, help='Where to store by product of the preprocessing',default='output')
parser.add_argument('-mo', '--outputMechanism', type=str, metavar='', required=False, help='Mechanism to write', default='qssa_final.inp')
parser.add_argument('-rv', '--forwardReactionsToRemove', type=str, metavar='', required=False, help='Detailed descriptions of forward reactions to remove', default='reac_forward_to_remove_detailed')
parser.add_argument('-r', '--forwardReactionsIDToRemove', type=str, metavar='', required=False, help='ID of forward reactions to remove', default='reac_forward_to_remove')
args = parser.parse_args()

# Set up filenames
mech_filename = args.inputMechanism
therm_filename = args.thermoFile
tran_filename = args.tranFile
forward_reac_remove_detailed_filename = args.forwardReactionsToRemove
forward_reac_remove_filename = args.forwardReactionsIDToRemove
non_qssa_species = args.nonQSSSpecies
outputFolder = args.outputFolder
output_mech_filename = args.outputMechanism



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
species_non_qssa = getListSpecies(non_qssa_species)
species_all = getListSpecies(app.inventory.mechanism)
species_qssa = list(set(species_all) - set(species_non_qssa))


# Write new mechanism
f = open(mech_filename,'r')
lines = f.readlines()
f.close()


# Get reactions to cancel
toCancel = identifyReactionsToCancel(mechanism,species_qssa)
print(toCancel)

outputFolder = outputFolder
mech_final_filename = output_mech_filename
f_mech = open(os.path.join(outputFolder,mech_final_filename),'w+')
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
                        break
                    else:
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
    f_mech.write(line_write)
f_mech.close()


#for ir, reaction in enumerate(mechanism.reaction()):
    


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



# Flag reactions to remove
f_forward_flag = open(os.path.join(outputFolder,forward_reac_remove_filename),'w+')
f_forward_detailed_flag = open(os.path.join(outputFolder,forward_reac_remove_detailed_filename),'w+')
line_iter = iter(lines)
for line in line_iter:
    line_write=line
    for ir in toCancel:
        if mechanism.reaction()[ir].equation().replace(' ','') in line:
            if 'f' in toCancel[ir] and 'r' in toCancel[ir]:
                break
            elif 'f' in toCancel[ir]:
                if mechanism.reaction()[ir].reversible:
                    if not mechanism.reaction()[ir].low and not mechanism.reaction()[ir].thirdBody:
                        f_forward_detailed_flag.write(str(mechanism.reaction()[ir]))
                        f_forward_detailed_flag.write('\n')
                        f_forward_flag.write(str(mechanism.reaction()[ir].id))
                        f_forward_flag.write('\n')
                        break
                    else:
                        f_forward_detailed_flag.write(str(mechanism.reaction()[ir]))
                        f_forward_detailed_flag.write('\n')
                        f_forward_flag.write(str(mechanism.reaction()[ir].id))
                        f_forward_flag.write('\n')
                        break
                else:
                    break
            elif 'r' in toCancel[ir]:
                    break
            else:
                break
f_forward_detailed_flag.close()
f_forward_flag.close()


if QSSCoupling(mechanism,species_qssa) == 0:
    print ('\n\n\nSUCCESS')
else:
    print ('\n\n\nFAILURE')

