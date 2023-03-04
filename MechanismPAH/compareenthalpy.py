

"""
This script evaluates the species content of the species: number of C, H, and O atoms
If any match between the mechanisms, L2 norm of the differences between the enthalpy values are checked
@author: Landon Owen
"""

import numpy as np
import matplotlib.pyplot as plt

Tlo = 300.
Thi = 3000.
# Enthalpy function
def h(T, coeff, Tmid):
    indx = 0
    if (T < Tmid):
        indx = 7
    curhval = coeff[indx+5] / T
    for j in range(5):
        curhval += coeff[indx+j] * T**j / float(j+1)
    return curhval
# Enthalpy function over array
def hv(T, coeff, Tmid):
    hval = np.zeros(len(T))
    for i in range(0,len(T)):
        Tv = T[i]
        hval[i] = h(Tv, coeff, Tmid)
    hval = np.array(hval)
    return hval
# Find number of C, H, and O in current species
def findspec(line, chstring):
    ccheck = chstring[0]
    hcheck = chstring[1]
    ocheck = chstring[2]
    numo = 1
    if (ocheck == "O   0"):
        numo = 0
    if (ccheck in line and hcheck in line):
        # Check for O
        if (ocheck in line):
            return True
        elif (numo == 0 and not "O" in line):
            return True
    return False

def getmidT(line):
    # Split at letter G
    lhalf = line.split("G")[1]
    # Split again based on spaces
    lvals = lhalf.split()
    if (len(lvals) == 4):
        return float(lvals[2])
    else:
        return 1000.
# It is important that both thermo files do not have unused species data in them
def compareenthalpy(pahthermo, newthermfile, plotH):
    files = [pahthermo, newthermfile]
    # These are the species most likely to have names that differ between files
    modspecs = ["C10H8", "C10H7", "C8H7", "C8H6", "C8H5", "C7H8", "C7H6O", "C7H7", "C6H5",  "C6H6", "C3H5"]
    namePAH = [] # Species name from PAH file
    nameMech = [] # Species name from mech file
    for curspec in modspecs:
        splitc = curspec.split("C")
        numc = int(splitc[1].split("H")[0])
        numh = int(splitc[1].split("H")[1][0])
        numo = 0
        if ("O" in curspec):
            numo = 1
        ocheck = "O   " + str(numo)
        ccheck = "C   " + str(numc)
        hcheck = "H   " + str(numh)
        if (numc > 9):
            ccheck = "C  " + str(numc)
        if (numh > 9):
            hcheck = "H  " + str(numh)
        chstring = [ccheck, hcheck, ocheck]
        numcols = [5, 5, 4]
        coeffs = []
        isfound = 0
        fspecname = []
        Tmid = []
        pahspec = ""
        pahcoeffs = []
        pahTmid = 0
        with open(pahthermo) as infile:
            for line in infile:
                spthere = findspec(line, chstring)
                if (spthere):
                    pahTmid = getmidT(line)
                    specname = line.split()[0]
                    pahspec = specname
                    for j in range(3):
                        line1 = next(infile)
                        for i in range(numcols[j]):
                            val1 = float(line1[0:15])
                            pahcoeffs.append(val1)
                            line1 = line1[15:]
        with open(newthermfile) as infile:
                for line in infile:
                    spthere = findspec(line, chstring)
                    if (spthere):
                        Tm = getmidT(line)
                        Tmid.append(Tm)
                        specname = line.split()[0]
                        isfound += 1
                        fspecname.append(specname)
                        coeffs.append([])
                        for j in range(3):
                            line1 = next(infile)
                            for i in range(numcols[j]):
                                val1 = float(line1[0:15])
                                coeffs[-1].append(val1)
                                line1 = line1[15:]
        # species is found in both files
        if (isfound >= 1):
            mechspecs = fspecname
            N = 200
            T = np.linspace(Tlo, Thi, N)
            cv = pahcoeffs
            hval1 = hv(T, cv, pahTmid)
            plt.plot(T, hval1)
            matches = 0
            for i in range(isfound):
                cv = coeffs[i]
                hval2 = hv(T, cv, Tmid[i])
                plt.plot(T, hval2)
                # Evaluate the L2 norm of the differences
                hdiff = abs(hval1 - hval2)
                l2norm = np.linalg.norm(hdiff, 2)
                if (l2norm > 1.):
                    print("WARNING: Enthalpy values for " + str(pahspec) + " and " + str(mechspecs[i]) + " differ significantly. We assume they are not the same species. However, further investigation is recommended. Try again with --plotH True.")
                else:
                    matches += 1
                    namePAH.append(pahspec)
                    nameMech.append(mechspecs[i])
            if (matches > 2):
                print("Found too many species for ")
            if (plotH):
                plt.show()
    return [namePAH, nameMech]
