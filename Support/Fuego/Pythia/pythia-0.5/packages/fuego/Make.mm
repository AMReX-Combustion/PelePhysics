# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = fuego

# directory structure

DIR_SRC = src
DIR_ETC = etc
DIR_MODULE = module
DIR_FUEGO = fuego
DIR_TESTS = tests
DIR_EXAMPLES = examples
DIR_APPLICATIONS = applications


BUILD_DIRS = \
    $(DIR_SRC) \
    $(DIR_MODULE) \
    $(DIR_FUEGO) \
    $(DIR_ETC) \
    $(DIR_APPLICATIONS) \

OTHER_DIRS = \
    $(DIR_TESTS) \
    $(DIR_EXAMPLES)

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

# targets

all: update

update: $(BUILD_DIRS)

release: tidy
	cvs release .

test: update
	(cd $(DIR_TESTS); $(MM) test)

.PHONY: $(DIR_SRC)
$(DIR_SRC):
	(cd $(DIR_SRC); $(MM))


.PHONY: $(DIR_ETC)
$(DIR_ETC):
	(cd $(DIR_ETC); $(MM))


.PHONY: $(DIR_MODULE)
$(DIR_MODULE):
	(cd $(DIR_MODULE); $(MM))


.PHONY: $(DIR_FUEGO)
$(DIR_FUEGO):
	(cd $(DIR_FUEGO); $(MM))


.PHONY: $(DIR_TESTS)
$(DIR_TESTS):
	(cd $(DIR_TESTS); $(MM))


.PHONY: $(DIR_EXAMPLES)
$(DIR_EXAMPLES):
	(cd $(DIR_EXAMPLES); $(MM))


.PHONY: $(DIR_APPLICATIONS)
$(DIR_APPLICATIONS):
	(cd $(DIR_APPLICATIONS); $(MM))


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
