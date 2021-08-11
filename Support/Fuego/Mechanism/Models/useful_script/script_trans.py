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

outlines = {}
for l in linesTR:
    l = l.strip()
    ts = l.split()
    if len(ts)>6 and ts[0] in spec_list:
        outlines[ts[0]] = '%-19s%1d%10.3f%10.3f%10.3f%10.3f%10.3f' % (ts[0],int(ts[1]),float(ts[2]),float(ts[3]),float(ts[4]),float(ts[5]),float(ts[6]))

# go through the list of species and find a match in the tran.dat file
# if several entries can be found the last one is kept
tran_list = []
for spec in spec_list:
    tran_list.append(outlines[spec])
    
# generate txt file with transport data, to append in the chem.inp
csv_file = 'TRANFILE_APPEND.txt'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)#,delimiter=' ',quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['TRANS ALL'])
    for kk in range(len(tran_list)):
        writer.writerow([tran_list[kk]])
    writer.writerow(['END'])

