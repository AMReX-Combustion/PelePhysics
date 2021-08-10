#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


def fastfind(fs, name):
    import re

    found = []
    worklist = [fs.root()]
    pattern = re.compile(name)

    for entry in worklist:
        for file in entry.children():
            if pattern.match(file.name):
                found.append(file)

        worklist += entry.subdirectories()

    return found


# version
__id__ = "$Id$"

#  End of file 
