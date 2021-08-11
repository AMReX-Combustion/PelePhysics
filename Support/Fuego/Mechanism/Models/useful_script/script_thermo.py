#!/usr/bin/env python

from __future__ import print_function
from builtins import range
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
    linesTH = infile1.readlines()
with open(file2, 'r') as infile2:
    linesMECH = infile2.readlines()

# find lines where species are defined
idx_beg = 1000000
idx_end = -1000000
for i, line in enumerate(linesMECH):
    if "SPECIES" in line:
        idx_beg = i
    if "END" in line and i > idx_beg:
        idx_end = i
        break

# generate a list of used species
spec_list = []
for i in range(idx_beg+1,idx_end):
    spec_list.extend(linesMECH[i].split())

print("Found ", len(spec_list), " species in mechanism")

# go through the list of species and find a match in the therm.dat file
# if several entries can be found the first one is kept
therm_list = []
for i, spec in enumerate(spec_list):
    therm_tmp = []
    print("Taking care of", spec.strip())
    for  j, line in enumerate(linesTH):
        if line.split():
            if line.split()[0] == spec.strip():
                #therm_tmp.append(line.strip())
                #therm_tmp.append(linesTH[j+1].strip())
                #therm_tmp.append(linesTH[j+2].strip())
                #therm_tmp.append(linesTH[j+3].strip())
                therm_tmp.append(line.split("\n")[0].split("\r")[0])
                therm_tmp.append(linesTH[j+1].split("\n")[0].split("\r")[0])
                therm_tmp.append(linesTH[j+2].split("\n")[0].split("\r")[0])
                therm_tmp.append(linesTH[j+3].split("\n")[0].split("\r")[0])
    if len(therm_tmp) > 4:
        print("  --> WARNING found additional entries in ",file1, ", taking the first one.")
    therm_list.append(therm_tmp[-4])
    therm_list.append(therm_tmp[-3])
    therm_list.append(therm_tmp[-2])
    therm_list.append(therm_tmp[-1])

# generate txt file with transport data, to append in the chem.inp
csv_file = 'THERMOFILE_APPEND.txt'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)#,delimiter=' ',quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['THERMO ALL'])
    for kk in range(len(therm_list)):
        writer.writerow([therm_list[kk]])
    writer.writerow(['END'])

