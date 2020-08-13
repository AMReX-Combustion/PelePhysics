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
ns = int(lines[14].split()[-1].split(':')[-1])
np = int(lines[15].split()[-1].split(':')[-1])

x_data = []
for i in range(86,nl-101,ns+1):
    x_data.append(float([x for x in lines[i].split()[-1].split('/') if x][-1]))

y_data = []
for j in range(0,ns+1):
    y_data.append([])
    for i in range(86+j,nl-101,ns+1):
        y_data[-1].append(float([x for x in lines[i].split()[-1].split('/') if x][-1]))


mech = sys.argv[2]
mech_file = os.environ['PELE_PHYSICS_HOME'] + '/Support/Fuego/Mechanism/Models/' + mech + '/mechanism.h'
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
lenx = len(x_data)
for i in range(lenx):
    f.write(' '.join([str(y_data[j][i]) for j in range(ns+1)])+'\n')
f.close()

print('pmf data writtent to ' + sys.argv[3])
