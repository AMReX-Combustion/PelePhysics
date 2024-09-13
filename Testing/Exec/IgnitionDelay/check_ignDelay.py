import sys

def readResult(filename):
    with open(filename,'r+') as f:
        line = f.readline()
        splitLine = line.split(',')
        temp = float(splitLine[0].split()[2])
        press = float(splitLine[1].split()[2])
        phi = float(splitLine[2].split()[2])
        tIgn = float(splitLine[3].split()[2])
        uncertainty = float(splitLine[3].split()[4])
    return tIgn, uncertainty



if __name__=="__main__":
    # Read ignition delay output
    tIgn,uncertainty = readResult('log')

    # Compare to theoretical value
    if abs(tIgn-0.0763)>15*uncertainty:
         sys.exit(1)
    else:
         sys.exit(0)
