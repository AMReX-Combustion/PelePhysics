# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = journal

#--------------------------------------------------------------------------
#

all: export

release: tidy
	cvs release .

update: clean
	cvs update .

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    Console.py \
    Device.py \
    Diagnostic.py \
    Entry.py \
    File.py \
    Index.py \
    IndexDebug.py \
    IndexError.py \
    IndexFirewall.py \
    IndexInfo.py \
    IndexWarning.py \
    Journal.py \
    NetRenderer.py \
    Renderer.py \
    TextFile.py \
    UDPDevice.py \
    __init__.py


export:: export-python-modules

# version
# $Id$

# End of file
