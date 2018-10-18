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
PACKAGE = cgi
PROJ_CLEAN += $(PROJ_CPPTESTS)

PROJ_PYTESTS = login.py
PROJ_TESTS = $(PROJ_PYTESTS)


#--------------------------------------------------------------------------
#

all: export

test:
	for test in $(PROJ_TESTS) ; do $${test}; done

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
    LoginForm.py \
    LoginPage.py \
    RegistrationForm.py \
    RegistrationPage.py \
    authorize.py \
    environment.py \
    login.py \
    register.py \
    settings.py \

export::
	-ssh ${WEB_HOST} mkdir -p ${WEB_DIR}
	-ssh ${WEB_HOST} $(RM) $(RMFLAGS) ${WEB_DIR}/*.pyc
	-ssh ${WEB_HOST} chmod a+rw ${WEB_DIR}/*
	scp -p ${EXPORT_WEB} ${WEB_HOST}:${WEB_DIR}


# version
# $Id$

# End of file
