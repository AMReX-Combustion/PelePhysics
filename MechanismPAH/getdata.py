

def addnewthermdata(pahthermofile, newthermofile, newPAHspec):
    """
    Add thermodynamic data for unrepresented species to the thermo file
    """
    newthermdata = []
    # First add thermodynamic data for new species
    with open(pahthermofile, 'r') as pfile:
        for line in pfile:
            spname = line.split()[0]
            if (spname in newPAHspec):
                newthermdata.append(line)
                for i in range(3):
                    newthermdata.append(next(pfile))
    with open(newthermofile, 'a+') as nfile:
        for newline in newthermdata:
            nfile.write(newline)
        nfile.write("END")
    return

def addnewtransdata(pahtransfile, newtransfile, newPAHspec):
    """
    Add transport data for unrepresented species to the thermo file
    """
    newtransdata = []
    with open(pahtransfile, 'r') as pfile:
        for line in pfile:
            spname = line.split()[0]
            if (spname in newPAHspec):
                newtransdata.append(line)
    with open(newtransfile, 'a+') as nfile:
        for newline in newtransdata:
            nfile.write(newline)
        nfile.write("END")
    return

def getaltreactnames(pahmechfile, modnamePAH, modnameMech):
    """
    This replaces species names in PAH reactions with the new names from the mech file
    """
    modpahreacts = []
    with open(pahmechfile, 'r') as pfile:
        for line in pfile:
            if (line[0:9] == "REACTIONS"):
                break
        for line in pfile:
            if (len(line.split()) > 0):
                if (line[0] != "!"):
                    outline = ""
                    # Remove the rate values
                    stoiline = line.split()[:-3]
                    modline = line
                    for n in range(len(stoiline)):
                        nv1 = modline.find(" ")
                        sval = modline[0:nv1]
                        if (sval in modnamePAH):
                            indx = modnamePAH.index(sval)
                            sval = modnameMech[indx]
                        outline += sval
                        modline = modline[nv1+1:]
                    if (not "END" in line):
                        modpahreacts.append(outline + modline)
    return modpahreacts

def updatemechfile(mechfile, newmechfile, pahmechfile, newPAHspec, modnamePAH, modnameMech):
    """
    Add new species in species list and remove transport data from mech file if present
    """
    fullmechfile = []
    with open(mechfile, 'r') as oldfile:
        for line in oldfile:
            if (line[0:7] == "SPECIES"):
                fullmechfile.append("SPECIES\n")
                newspecstr = ""
                for newspec in newPAHspec:
                    newspecstr += newspec + " "
                newspecstr += "\n"
                fullmechfile.append(newspecstr)
            elif (line[0:9] == "TRANS ALL"):
                # Ignore transport properties
                while (not "END" in line):
                    line = next(oldfile)
            elif (len(line) > 0):
                fullmechfile.append(line)
    # Remove the last END
    fullmechfile = fullmechfile[:-1]
    # Add line signifying the start of the PAH mechanism
    fullmechfile.append("! Start of PAH module\n")
    # Get reaction PAH reaction data, rename species where necessary
    # Note: This call is specific to the PAH mech file and not all mech files
    modpahreacts = getaltreactnames(pahmechfile, modnamePAH, modnameMech)
    with open(newmechfile, 'w') as mfile:
        for line in fullmechfile:
            mfile.write(line)
        for line in modpahreacts:
            mfile.write(line)
        mfile.write("END\n")
    return [modpahreacts, fullmechfile]


