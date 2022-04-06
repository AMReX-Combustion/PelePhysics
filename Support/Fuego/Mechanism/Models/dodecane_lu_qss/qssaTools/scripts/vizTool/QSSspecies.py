def getListSpecies(filename):
    f = open(filename,'r+')
    lines = f.readlines()
    f.close()
    # Get species names    
    species = []
    startRead = False
    for line in lines:
        if line.startswith('SPECIES'):
           startRead = True
        if line.startswith('END') and startRead:
           break
        if not line.startswith('SPECIES') and startRead:
           species += line.split()
    return species
