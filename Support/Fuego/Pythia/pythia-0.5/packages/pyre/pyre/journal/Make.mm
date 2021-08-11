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
PACKAGE = journal


#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    Channel.py \
    ChannelFacility.py \
    Console.py \
    Daemon.py \
    Device.py \
    DeviceFacility.py \
    File.py \
    Journal.py \
    JournalFacility.py \
    NetRenderer.py \
    Remote.py \
    Renderer.py \
    RendererFacility.py \
    Server.py \
    __init__.py


export:: export-package-python-modules

# version
# $Id$

# End of file
