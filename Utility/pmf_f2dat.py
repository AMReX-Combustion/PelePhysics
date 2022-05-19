import sys
import os

if len(sys.argv) != 4:
    print(len(sys.argv))
    sys.exit('Usage: python ' + sys.argv[0] + ' <input pmf.f filename> <mechID from PP> <out pmf.dat filename>')

print('reading pmf fortran file ' + sys.argv[1])

try:
    f = open(sys.argv[1])
except Exception as ex:
    sys.exit(ex)
    
lines = f.readlines()
f.close()

nl = len(lines)

for l in lines:
    split = l.strip().split(':')
    if 'umber of states' in split[0]:
        ns = int(split[-1])
        break

for l in lines:
    split = l.strip().split(':')
    if 'umber of x data points' in split[0]:
        np = int(split[-1])
        break

print('Number of states per point',ns)
print('Number of points',np)

in_dat = False
nl_hdr = 0
for i in range(nl):
    s = lines[i].split()
    is_nondat = len(s)==1 or not (len(s) > 2 and s[0] == 'data' and s[1] == 'x_data(1)')
    if not in_dat:
        if is_nondat:
            nl_hdr = nl_hdr+1
        else:
            in_dat = True

x_data = []
for i in range(nl_hdr,nl+1,ns+1):
    if lines[i].split()[0] != 'data':
        break
    x_data.append(float([x for x in lines[i].split()[-1].split('/') if x][-1]))
np = len(x_data)

y_data = []
for j in range(1,ns+1):
    y_data.append([])
    for i in range(nl_hdr+j,nl+1,ns+1):
        if lines[i].split()[0] != 'data':
            break
        y_data[-1].append(float([x for x in lines[i].split()[-1].split('/') if x][-1]))
    if len(y_data[-1]) != np:
        sys.exit('Different number of values read for y_data than x_data')

mech = sys.argv[2]
mech_file = os.environ['PELE_PHYSICS_HOME'] + '/Support/Mechanism/Models/' + mech + '/mechanism.H'
print('grabbing species list from ' + mech_file)
try:
    f = open(mech_file)
except Exception as ex:
    sys.exit(ex)

lines = f.readlines()
f.close()

species = []
for l in lines:
    tok = l.split()
    if len(tok)==3 and '_ID' in tok[1]:
        species.append('"' + tok[1][:-3] + '"')
print('Found '+str(len(species))+' species names to use')
if len(species) + 3 != ns:
    sys.exit('Number of species not compatible with fortran file')

try:
    f = open(sys.argv[3],'w')
except Exception as ex:
    sys.exit(ex)

f.write('VARIABLES = "X" "temp" "u" "rho" ' + ' '.join(species)+'\n')
f.write('ZONE I=' + str(np) + ' FORMAT=POINT\n')
for i in range(np):
    f.write(str(x_data[i])+' '+' '.join([str(y_data[j][i]) for j in range(ns)])+'\n')
f.close()

print('pmf data writtent to ' + sys.argv[3])
