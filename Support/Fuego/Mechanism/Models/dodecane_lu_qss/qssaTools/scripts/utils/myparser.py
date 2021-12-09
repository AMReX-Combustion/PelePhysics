import sys

def parseInputFile():

    try:
        input_fileName = sys.argv[1]
    except IndexError:
        print("No input file name found, assuming it is 'input'")
        input_fileName = 'input'

    # ~~~~ Parse input
    inpt = {}
    f = open( input_fileName )
    data = f.readlines()
    for line in data:
        if ':' in line:
            key, value = line.split(":")
            inpt[key.strip()] = value.strip()
    f.close()
   
    return inpt
