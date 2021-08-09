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

import os
import sys

GCC_DIR = "/home/tools/gcc-3.1"
PYTHIA_DIR = "/home/tools/pythia-0.1"

def configure():
    fixIO()
    fixEnvironment()
    installPackages()
    return


def fixIO():
    sys.stderr = sys.stdout
    
    import cgitb
    cgitb.enable()
    
    return


def fixEnvironment():
    # LD_LIBRARY_PATH
    LD_LIBRARY_PATH = os.environ.get("LD_LIBRARY_PATH", "")
    if LD_LIBRARY_PATH:
        LD_LIBRARY_PATH += ":"

    LD_LIBRARY_PATH += PYTHIA_DIR + "/lib" + ":" + GCC_DIR + "/lib"
    os.environ["LD_LIBRARY_PATH"] = LD_LIBRARY_PATH
    
    return


def installPackages():
    sys.path = [PYTHIA_DIR] + sys.path

    # installJournal()

    return


def installJournal():
    import journal
    return


# version
__id__ = "$Id$"

#  End of file 
