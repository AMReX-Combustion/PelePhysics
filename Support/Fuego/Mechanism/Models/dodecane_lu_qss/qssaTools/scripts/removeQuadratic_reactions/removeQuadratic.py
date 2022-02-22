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


# Get reactions to cancel
toCancel = identifyReactionsToCancel(mechanism,species_qssa)
print(toCancel)

# Write new mechanism
f = open(mech_filename,'r')
lines = f.readlines()
f.close()

outputFolder = inpt['outputFolder'] 
mech_final_filename = inpt['mech_final']
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
            #elif 'f' in toCancel[ir]:
            #    if mechanism.reaction()[ir].reversible:
            #        if not mechanism.reaction()[ir].low and not mechanism.reaction()[ir].thirdBody:
            #            decompline = line.split()
            #            print(decompline)
            #            reac_trans_prod = decompline[0].split('<=>')
            #            line_write = reac_trans_prod[1] + '=>' + reac_trans_prod[0]
            #            for entry in decompline[1:]:
            #                line_write += '\t' + entry
            #            line_write += '\n'
            #            print(mechanism.reaction()[ir].equation())
            #            print(line_write)
            #            break
            #        else:
            #            line_write=''
            #            next(line_iter)
            #            next(line_iter)
            #            next(line_iter)
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


if QSSCoupling(mechanism,species_qssa) == 0:
    print ('\n\n\nSUCCESS')
else:
    print ('\n\n\nFAILURE')

