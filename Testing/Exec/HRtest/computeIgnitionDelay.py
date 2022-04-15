import sys
sys.path.append('utils')
import argparse
import fileio as io
import numpy as np


def parseInputFile(input_filename):
    # ~~~~ Parse input
    inpt = {}
    f = open( input_filename )
    data = f.readlines()
    for line in data:
        if '=' in line:
            splitLine = line.split("=")
            if len(splitLine)==2:
                key, value = splitLine
            else:
                key, value = splitLine[:2]
            inpt[key.strip()] = value.strip()
    f.close()

    return inpt

def writeNewInput(oldInput,newInput,dt,ndt):
    f = open(oldInput,'r+')
    lines = f.readlines()
    f.close()

    f = open(newInput,'w+')
    for line in lines:
        if line.startswith('ode.dt'):
            f.write('ode.dt = %.10f' % (dt))
            f.write('\n')
        elif line.startswith('ode.ndt'):
            f.write('ode.ndt = %d ' % (ndt))
            f.write('\n')
        else:
            f.write(line)
    f.close()

    return

# CLI
parser = argparse.ArgumentParser(description='Compute ignition delay')
group = parser.add_mutually_exclusive_group()
parser.add_argument('-f', '--inputFile', type=str, metavar='', required=False, help='Pele input file', default='')
group.add_argument('-est','--estimate', action='store_true', help='first estimate')
group.add_argument('-ref','--refine', action='store_true', help='Second estimate')
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-v','--verbose', action='store_true', help='verbose')
group2.add_argument('-q','--quiet', action='store_true', help='quiet')
args = parser.parse_args()


# Read Pele input file
inpt = parseInputFile(args.inputFile)
dt = float(inpt['ode.dt'])
ndt = float(inpt['ode.ndt'])
press = float(inpt['hr.press'])
temp = float(inpt['hr.t0'])
phi = float(inpt['hr.equiv_ratio'])


# Process result
filePP = 'PPreaction.txt'
A1 = io.readMultiColFile(filePP)
ignitionDelayPP = (dt/ndt)*np.amin(np.argwhere(A1[:,1]>(A1[0,1]+400)))


if args.estimate:
    if args.verbose: 
        print("Rough estimate: T = %.3g K , P = %.3g atm, Phi = %.3g,  tIgn = %.3g s " % (temp,press/101325.0,phi,ignitionDelayPP))
    # Write new input file
    writeNewInput(args.inputFile,'inputs/inputs.0d_refine',ignitionDelayPP*1.1,ndt*2)

if args.refine:
    print("T = %.3g K , P = %.3g atm, Phi = %.3g,  tIgn = %.3g s " % (temp,press/101325.0,phi,ignitionDelayPP))
