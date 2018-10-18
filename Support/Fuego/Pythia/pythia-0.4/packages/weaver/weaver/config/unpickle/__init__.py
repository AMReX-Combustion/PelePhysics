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


def readUserConfiguration():
    import os

    try:
        filename = os.path.join(os.environ["HOME"], ".weaver")
    except KeyError:
        return None
    
    try:
        userfile = file(filename, "r")
    except IOError:
        import journal
        info = journal.info("weaver.config")
        info.log("could not read '%s'" % filename)
        return None

    options = readConfiguration(userfile)
    return options


def readConfiguration(file):

    import pyre.xml
    from Configuration import Configuration

    p = pyre.xml.parser()
    p.parse(file, Configuration(file.name))

    options = p.document().options()

    return options

# version
__id__ = "$Id$"

#  End of file 
