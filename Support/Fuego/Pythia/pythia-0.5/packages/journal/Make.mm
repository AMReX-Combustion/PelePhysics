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

PROJECT = journal

# directory structure

DIR_LIB = lib
DIR_MODULE = module
DIR_JOURNAL = journal
DIR_TESTS = tests
DIR_EXAMPLES = examples


BUILD_DIRS = \
    $(DIR_LIB) \
    $(DIR_MODULE) \
    $(DIR_JOURNAL) \

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

.PHONY: $(DIR_LIB)
$(DIR_LIB):
	(cd $(DIR_LIB); $(MM))


.PHONY: $(DIR_MODULE)
$(DIR_MODULE):
	(cd $(DIR_MODULE); $(MM))


.PHONY: $(DIR_JOURNAL)
$(DIR_JOURNAL):
	(cd $(DIR_JOURNAL); $(MM))


.PHONY: $(DIR_TESTS)
$(DIR_TESTS):
	(cd $(DIR_TESTS); $(MM))


.PHONY: $(DIR_EXAMPLES)
$(DIR_EXAMPLES):
	(cd $(DIR_EXAMPLES); $(MM))


tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse
	$(RM) $(RMFLAGS) $(BLD_TMPDIR)/$(PROJECT) $(BLD_LIBDIR)/$(PROJECT)

distclean::
	BLD_ACTION="distclean" $(MM) recurse


# version
# $Id$

#
# End of file
