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

PROJECT = weaver

# directory structure

DIR_APPLICATIONS = applications
DIR_WEAVER = weaver
DIR_TESTS = tests
DIR_EXAMPLES = examples
DIR_HTML = html


BUILD_DIRS = \
    $(DIR_WEAVER) \
    $(DIR_APPLICATIONS) \

OTHER_DIRS = \
    $(DIR_HTML) \
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


.PHONY: $(DIR_APPLICATIONS)
$(DIR_APPLICATIONS):
	(cd $(DIR_APPLICATIONS); $(MM))


.PHONY: $(DIR_WEAVER)
$(DIR_WEAVER):
	(cd $(DIR_WEAVER); $(MM))


.PHONY: $(DIR_TESTS)
$(DIR_TESTS):
	(cd $(DIR_TESTS); $(MM))


.PHONY: $(DIR_EXAMPLES)
$(DIR_EXAMPLES):
	(cd $(DIR_EXAMPLES); $(MM))


.PHONY: $(DIR_HTML)
$(DIR_HTML):
	(cd $(DIR_HTML); $(MM))


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
