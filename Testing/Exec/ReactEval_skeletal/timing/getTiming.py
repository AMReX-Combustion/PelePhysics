import numpy as np
import sys

def getTiming(logFile):
    f=open(logFile,'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if " >> React::main()" in line:
             time=float(line.split()[-1])        
             return time

def writeTiming(time,filename):
    f=open(filename,'a+')
    f.write(str(time)+'\n')
    return  

logFile = sys.argv[1]
timeFile = sys.argv[2]


time = getTiming(logFile)
writeTiming(time,timeFile)

