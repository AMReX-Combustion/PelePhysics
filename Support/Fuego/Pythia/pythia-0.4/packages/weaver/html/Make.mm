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

PROJECT = weaver
PACKAGE = html

#--------------------------------------------------------------------------
#

DIR_IMAGES = images

OTHER_DIRS = $(DIR_IMAGES)
RECURSE_DIRS = $(OTHER_DIRS)

all: export

release: tidy
	cvs release .

update: clean
	cvs update .

#--------------------------------------------------------------------------
#
# export

WEB_HOST = nectar.cacr.caltech.edu
WEB_DIR = /home/web/${PROJECT}/${PACKAGE}
EXPORT_WEB = \
    index.html \
    test-style.html \


export::
	-ssh ${WEB_HOST} mkdir -p ${WEB_DIR}
	-ssh ${WEB_HOST} chmod a+rw ${WEB_DIR}/*
	scp -p ${EXPORT_WEB} ${WEB_HOST}:${WEB_DIR}


#--------------------------------------------------------------------------
#
# clean up
tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse
	$(RM) $(RMFLAGS) $(BLD_TMPDIR)/$(PROJECT) $(BLD_LIBDIR)/$(PROJECT)

distclean::
	BLD_ACTION="distclean" $(MM) recurse


# version
# $Id$

# End of file
