from FMC import FMC
import fuego
from QSSspecies import getListSpecies 
from coupling import *
import sys
sys.path.append('scripts/utils')
import myparser
import os



# Parse input
inpt = myparser.parseInputFile()

# Set up filenames
mech_filename = inpt['mech']
therm_filename = inpt['therm']
tran_filename = inpt['tran']
forward_reac_remove_detailed_filename = inpt['reac_forward_remove_detailed']
forward_reac_remove_filename = inpt['reac_forward_remove']

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


# Write new mechanism
f = open(mech_filename,'r')
lines = f.readlines()
f.close()


# Get reactions to cancel
toCancel = identifyReactionsToCancel(mechanism,species_qssa)
print(toCancel)

outputFolder = inpt['outputFolder'] 
mech_final_filename = inpt['mech_final']
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

