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
    # Two arguments: intended correct value and uncertainty factor
    if len(sys.argv) > 1:
        correct_value = float(sys.argv[1])
        uncertainty_factor = float(sys.argv[2])
    else:
        correct_value = 0.0763
        uncertainty_factor = 10.0

    # Read ignition delay output
    tIgn,uncertainty = readResult('log')

    # Compare to theoretical value
    if abs(tIgn-correct_value)>uncertainty_factor*uncertainty:
         sys.exit(1)
    else:
         sys.exit(0)
