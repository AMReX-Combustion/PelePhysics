import argparse
import numpy as np
import sys

def parseInputFile(input_filename):

    inpt = {}
    with open(input_filename,'r+') as f:
        data = f.readlines()
        for line in data:
            # Remove comments from the line
            commentStart = line.find('#')
            if commentStart>-1:
                line = line[:commentStart]
            # Store parameters in dict
            if '=' in line:
                splitLine = line.split("=")
                if len(splitLine)==2:
                    key, value = splitLine
                else:
                    key, value = splitLine[:2]
                inpt[key.strip()] = value.strip()

    return inpt

def writeNewInput(oldInput,newInput,dt,ndt):
    
    with open(oldInput,'r+') as f: 
        lines = f.readlines()

    with open(newInput,'w+') as f: 
        for line in lines:
            if line.startswith('ode.dt'):
                f.write('ode.dt = %.10f' % (dt))
                f.write('\n')
            elif line.startswith('ode.ndt'):
                f.write('ode.ndt = %d ' % (ndt))
                f.write('\n')
            else:
                f.write(line)


def readMultiColFile(filename,headerSize=0):
    with open(filename,'r') as f:
        lines=f.readlines()
    # Figure out number of column
    A=lines[headerSize].split()
    # Allocate array
    Array = np.zeros((len(lines)-headerSize,len(A)))
    counter = 0
    for line in lines[headerSize:]:
        A=line.split()
        for i in range(len(A)):
            Array[counter,i] = float(A[i])
        counter=counter+1
    return Array




if __name__=="__main__":

    # CLI
    parser = argparse.ArgumentParser(description='Compute ignition delay')
    group = parser.add_mutually_exclusive_group()
    parser.add_argument('-f', '--inputFile', type=str, metavar='', required=False, help='Pele input file', default='')
    parser.add_argument('-t','--test', action='store_true', help='testing mode, only for Github test')
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
    A1 = readMultiColFile(filePP)
    ignitionDelayPP = (dt/ndt)*np.amin(np.argwhere(A1[:,1]>(A1[0,1]+400)))
    
    
    if args.estimate:
        if args.verbose: 
            logMessage = f"Rough estimate: T = {temp:.3g} K , P = {press/101325.0:.3g} atm, Phi = {phi:.3g},  tIgn = {ignitionDelayPP:.3g} +/- {dt/ndt:.3g} s"
            print(logMessage)
        # Write new input file
        writeNewInput(args.inputFile,'inputs/inputs.0d_refine',ignitionDelayPP*1.1,ndt*2)
    
    if args.refine: 
        logMessage = f"T = {temp:.3g} K , P = {press/101325.0:.3g} atm, Phi = {phi:.3g},  tIgn = {ignitionDelayPP:.3g} +/- {dt/ndt:.3g} s"
        print(logMessage)
        print("Equivalence ratio calculation assumed you are using Ndodecane")
        print("Change GPU_misc.H if you are not")
        if args.test:
            if abs(ignitionDelayPP-0.0763)>10*dt/ndt:
                 sys.exit(1)
            else:
                 sys.exit(0)
