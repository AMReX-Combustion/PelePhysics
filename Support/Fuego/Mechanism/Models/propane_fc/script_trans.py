import numpy as np
import sys
import csv

####################################################################
#
#  USAGE: "script_trans.py chem.inp tran.dat"
#         It will return a file labelled TRANFILE_APPEND.txt. 
#         The user should copy all the lines in the chem.inp file,
#         right above the reactions definitions for example.
#
####################################################################

# open files
file1 = sys.argv[2]
file2 = sys.argv[1]
with open(file1, 'r') as infile1:
    linesTR = infile1.readlines()
with open(file2, 'r') as infile2:
    linesGRI = infile2.readlines()

# find lines where species are defined
idx_beg = 1000000
idx_end = -1000000
for i, line in enumerate(linesGRI):
    if "SPECIES" in line:
        idx_beg = i
    if "END" in line and i > idx_beg:
        idx_end = i
        break

# generate a list of used species
spec_list = []
for i in range(idx_beg+1,idx_end):
    spec_list.extend(linesGRI[i].split())

print("Found ", len(spec_list), " species in mechanism")

# go through the list of species and find a match in the tran.dat file
# if several entries can be found the last one is kept
tran_list = []
for i, spec in enumerate(spec_list):
    tran_tmp = []
    print("Taking care of", spec.strip())
    for  j, line in enumerate(linesTR):
        if line.split():
            if line.split()[0] == spec.strip():
                #print("Found it in tranfile", (line.split("!")[0]).strip())
                tran_tmp.append((line.split("!")[0]).strip()+"          ")
                #tran_list.append((line.split("!")[0]).strip()+"          ")
    if len(tran_tmp) > 1:
        print("  --> WARNING found additional entries in ",file1, ", taking the last one.")
    tran_list.append(tran_tmp[-1])

#if len(tran_list)!=len(spec_list):
#    print("PROBLEM !! Len tran_list is ",len(tran_list), " when len spec_list is ",len(spec_list))
#    STOP

# generate txt file with transport data, to append in the chem.inp
csv_file = 'TRANFILE_APPEND.txt'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)#,delimiter=' ',quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['TRANS ALL'])
    for kk in range(len(tran_list)):
        writer.writerow([tran_list[kk]])
    writer.writerow(['END'])

