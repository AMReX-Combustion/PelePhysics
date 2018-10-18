# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2001  All Rights Reserved
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def

PROJECT = cv-driver
PROJ_LIB = $(BLD_LIBDIR)/$(PROJECT).$(EXT_LIB)
PROJ_CLEAN += *~ core cv

MECHANISM = $(PYTHIA_DIR)/share/mechanisms/HydrogenOxygen.ck2

# Project source files
PROJ_SRCS = \
    ddebdf.$(EXT_F77) \
    burn_cv.$(EXT_F77) \
    fuego.$(EXT_C)


LIBS = $(PROJ_LIB) $(EXTERNAL_LIBS) -lm


all: $(PROJ_LIB) cv


cv: driver.$(EXT_C)
	$(CC) -o cv driver.$(EXT_C) $(CFLAGS) $(LCFLAGS) $(LIBS)


fuego.$(EXT_C): fuego.py $(MECHANISM)
	python fuego.py --file=fuego.$(EXT_C) --mechanism=$(MECHANISM) --output=c


# version
# $Id$

#
# End of file
