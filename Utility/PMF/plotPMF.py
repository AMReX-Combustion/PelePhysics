"""Plot data from PMF files."""

import sys

import matplotlib.pyplot as plt
import numpy as np

try:
    pmf = np.genfromtxt(sys.argv[1], skip_header=2)
except IOError:
    print("There was an error reading data from the file: " + sys.argv[1])
    sys.exit()

try:
    f = open(sys.argv[1])
except IOError:
    print("There was an error open the file: " + sys.argv[1])
    sys.exit()

vars = f.readline().strip().replace('"', "").replace(",", "").split()[2:]
if len(vars) != pmf.shape[1]:
    print("Variable list not consistent with data in " + sys.argv[1])

var_map = {}
for i, var in enumerate(vars):
    var_map[var] = i

print(var_map)

vid = var_map["temp"]
if len(sys.argv) > 2:
    vid = int(sys.argv[2])
plt.plot(pmf[:, 0], pmf[:, vid], linewidth=1.0, color="black")
plt.xlabel("Position (cm)", fontsize=14)
plt.ylabel(vars[vid], fontsize=14)
plt.grid(True)
plt.show()
