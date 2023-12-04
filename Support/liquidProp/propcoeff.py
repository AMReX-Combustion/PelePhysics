#!/usr/bin/env python

# This code reads in an NIST property file for a
# condensed or saturated phase of a species
# It fits the data to data depending on the parameter
# Should have files for rho, mu, and lambda (thermal conductivity)
# They should be in a single directory and called rho.dat, mu.dat, and lambda.dat
# Currently only works for hydrocarbons

# Usage:
# python propfit.py --species NC10H22 --file_loc decane --units MKS --vars rho

import sys
import csv
import os
import numpy.polynomial.polynomial as nppp
import numpy as np
from scipy.optimize import curve_fit
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("--species", help="Species name", metavar="NC10H22", type=str,required=True)
parser.add_argument("--file_loc", help="Location of data files. Files should be called rho.dat, mu.dat, and/or lambda.dat", type=str,default="None")
parser.add_argument("--units", help="Units, either MKS or CGS", type=str,default="MKS")
parser.add_argument("--vars", help="Which variables to fit, ex. mu lambda rho", default="None",nargs='+', type=str)
arg_string = sys.argv[1:]
args = parser.parse_args()
Tmin = 280.

def mufit(T, a, b, c, d):
    return a + ((d / T + c) / T + b) / T

def genfit(T, a, b, c, d):
    return a + T * (b + T * (c + d * T))
def varfitdata(infile, func, conv):
    T = []
    prop = []
    if (not os.path.exists(infile)):
        error = "file " + infile + " does not exist"
        raise ValueError(error)
    origT = []
    origProp = []
    with open(infile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        k = 0
        for row in csv_reader:
            if (k > 2):
                if (float(row[0]) > Tmin):
                    origT.append(float(row[0]))
                    origProp.append(float(row[1]) * conv)
                    # Check if uncertainty is reasonable
                    ucq = float(row[-1])
                    kv = float(row[1])
                    if (ucq/kv < 0.1):
                        T.append(float(row[0]))
                        prop.append(float(row[1]) * conv)
            k += 1
    # If all points have high error, just use all of them
    if (len(T) < 4):
        T = np.array(origT)
        prop = np.array(origProp)
    else:
        T = np.array(T)
        prop = np.array(prop)
    coeff, res = curve_fit(func, T, prop)
    vfit = func(T, *coeff)
    aprop = prop
    return [coeff, T, aprop, vfit]

def psat(T, a, b, c):
    expn = a - b / (T + c)
    return expn
def fitpsat(infile):
    T = []
    prop = []
    if (not os.path.exists(infile)):
        error = "file " + infile + " does not exist"
        raise ValueError(error)
    with open(infile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        k = 0
        for row in csv_reader:
            if (k > 2):
                # Check if uncertainty is reasonable
                ucq = float(row[-1])
                kv = float(row[1])
                if (float(row[0]) > Tmin):# and ucq/kv < 0.1):
                    T.append(float(row[0]))
                    prop.append(np.log10(float(row[1])*0.01))
            k += 1
    T = np.array(T)
    prop = np.array(prop)
    coeff, res = curve_fit(psat, T, prop)
    Tv = np.linspace(T[0], T[-1],len(T))
    vfit = 10.**(psat(Tv, *coeff))
    aprop = 10.**prop
    return [coeff, T, aprop, vfit]

def mufit(T, a, b, c, d):
    return a + ((d / T + c) / T + b) / T

if (len(arg_string) == 0):
    parser.print_help()
    sys.exit()

defvarnames = ["lambda", "mu", "rho", "cp", "sigma"]
varnames = []

if (args.vars == "None"):
    for i in range(len(defvarnames)):
        varnames.append(defvarnames[i])
else:
    for i in range(len(args.vars)):
        varnames.append(args.vars[i])
numvar = len(varnames)
speciesname = args.species
# Extract species to get molecular weight
# Currently unused but could be used if C_p or enthalpy comes in
splitc = speciesname.split("C")
[numcst, numhst] = splitc[1].split("H")
numc = int(numcst)
numh = int(numhst)
mw = numc * 12. + numh
# convert to mks of kg/mol
mw /= 1000.
if (args.file_loc == "None"):
    args.file_loc = speciesname

# For converting to CGS
mu_conv = 1.E5 / (100.**2)
lambda_conv = 1.E7 / 1000.
rho_conv = 1000. / (100.**3)
p_conv = 10.
s_conv = 1000.
cp_conv = 1.E7 / 1000.
cgs_conv = []

# Psat and mu are different
for i in range(numvar):
    if (varnames[i] == "mu"):
        cgs_conv.append(mu_conv)
    elif (varnames[i] == "lambda"):
        cgs_conv.append(lambda_conv)
    elif (varnames[i] == "rho"):
        cgs_conv.append(rho_conv)
    elif (varnames[i] == "psat"):
        cgs_conv.append(p_conv)
    elif (varnames[i] == "cp"):
        cgs_conv.append(cp_conv)
    elif (varnames[i] == "sigma"):
        cgs_conv.append(s_conv)
    else:
        errorstatement = "Variable "+ varnames[i] + " not recognized"
        raise ValueError(errorstatement)

coeffs = []
for f in range(len(varnames)):
    infile = args.file_loc + "/" + varnames[f] + ".dat"
    if (varnames[f] == "psat"):
        [p, T, prop, vfit] = fitpsat(infile)
        if (args.units == "CGS"):
            p = np.append(p, 1.E6)
        else:
            p = np.append(p, 1.E5)
        plt.plot(T, vfit)
        plt.scatter(T, prop)
        plt.show()
    else:
        conv = 1.
        if (varnames[f] == "mu"):
            func = mufit
        else:
            if (varnames[f] == "cp"):
                conv = 1. / mw
            func = genfit
        [p, T, prop, vfit] = varfitdata(infile, func, conv)
        if (args.units == "CGS"):
            for i in range(len(p)):
                p[i] *= cgs_conv[f]
        plt.plot(T, vfit)
        plt.scatter(T, prop)
        plt.show()
    coeffs.append([])
    for i in range(len(p)):
        coeffs[f].append(p[i])
for f in range(len(varnames)):
    curcoeffs = ""
    for i in range(len(coeffs[f])):
        curcoeffs += "{:g} ".format(coeffs[f][i])
    print("particles." + args.file_loc + "_" + varnames[f] + " = " + curcoeffs)
