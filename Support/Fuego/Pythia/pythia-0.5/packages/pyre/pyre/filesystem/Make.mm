# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = pyre
PACKAGE = filesystem


#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    BlockDevice.py \
    CharacterDevice.py \
    Directory.py \
    Entry.py \
    Explorer.py \
    fastfind.py \
    File.py \
    FileSystem.py \
    Finder.py \
    Inspector.py \
    Link.py \
    NamedPipe.py \
    SimpleRenderer.py \
    Socket.py \
    TreeRenderer.py \
    __init__.py


export:: export-package-python-modules

# version
# $Id$

# End of file
