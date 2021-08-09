#!/usr/bin/env python

import sys
import os

####################################################################
#
#  USAGE: "python script_generate_mechs.py listOfMechs"
#         This script will go into each folder as listed in
#         listOfMechs and launch the shell script that regenerates
#         a cpp mechanism with Fuego
#  NOTE: You must source the file located at $PELE_PHYSICS/Support/Fuego/Pythia/setup.sh
#        Before launching this script
#
####################################################################

ListOfFolders = sys.argv[1]
try:
    # open list of mechs to update
    with open(ListOfFolders, 'r') as infile:
        Folders = infile.readlines() 
except:
    # otherwise deal with just one mech
    Folders = [sys.argv[1]]

# Root of Models
root_dir = os.getcwd()

# Warning blabla
print(" --> PLEASE READ !! ") 
print("     -----------    ")
print("You should NOT run this script unless you understand what it does. ")
print("The usage is 'python script_generate_mechs.py listOfMechs' where the listOfMechs is a simple list of mechanisms ")
print("that needs to be updated. ")
print("Outputs for EACH mechanisms are redirected in the relevant folder and should be ")
print("carefully examined to understand errors and fix them. ")
print("A badly translated mechanism will lead to trash simulations.\n")
print("You have been warned. Type 'yes' if you understand and want to pursue. Anything else will abort: ")
areYouSure = raw_input()

if (areYouSure == "yes"):
    # make cpp for requested mechs
    for i, mech in enumerate(Folders):
        os.chdir(root_dir+"/"+mech.strip()) 
        curr_dir = os.getcwd()
        print("#---------#---------#---------#---------#---------#") 
        print("#     Taking care of " +mech.strip()+ "...        #")
        print("#---------#---------#---------#---------#---------#") 
        os.system("./make-mechanism.sh > out_"+mech.strip())

else:
    # abort
    print("aborting")
