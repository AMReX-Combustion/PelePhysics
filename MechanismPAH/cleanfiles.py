
import sys

# This routine reads in the thermo and transport files and removes any unused species

def cleanfiles(speciesnames, thermfile, transfile, newthermfile, newtransfile):
    indatfile = thermfile
    outdatfile = newthermfile
    foundspecs = 0
    foundnames = []
    with open(indatfile) as infile:
        with open(outdatfile, 'w') as outfile:
            # Copy first two lines over
            line = infile.readline()
            outfile.write(line)
            line = infile.readline()
            outfile.write(line)
            for line in infile:
                if (len(line.split()) > 0):
                    if (line[0] != '!'):
                        spname = line.split()[0]
                        for sp in speciesnames:
                            if (sp == spname):
                                foundnames.append(sp)
                                foundspecs += 1
                                outfile.write(line)
                                for i in range(3):
                                    line = next(infile)
                                    outfile.write(line)
    if (foundspecs != len(speciesnames)):
        for sp in speciesnames:
            if (not sp in foundnames):
                print(sp + " not found in thermo file")
        print("Not all species were found in thermo file")
        sys.exit()
    print("found " + str(foundspecs) + " out of " + str(len(speciesnames)) + " for thermo data")
    indatfile = transfile
    outdatfile = newtransfile
    foundspecs = 0
    foundnames = []
    with open(indatfile) as infile:
        with open(outdatfile, 'w') as outfile:
            outfile.write("TRANS ALL\n")
            for line in infile:
                if (len(line.split()) > 0):
                    if (line[0] != '!'):
                        spname = line.split()[0]
                        for sp in speciesnames:
                            if (sp == spname):
                                foundspecs += 1
                                outfile.write(line)
    if (foundspecs != len(speciesnames)):
        for sp in speciesnames:
            if (not sp in foundnames):
                print(sp + " not found in transport file")
        print("Not all species were found in transport file")
        sys.exit()
    print("Found " + str(foundspecs) + " out of " + str(len(speciesnames)) + " for transport data")
    print("Successfully made clean versions of thermodynamic and transport files")
