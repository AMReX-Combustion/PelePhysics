# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = pyre
PROJ_TIDY += $(PROJ_CPPTESTS) $(PROJ_CTESTS)
PROJ_CLEAN = $(PROJ_TIDY)

PROJ_PYTESTS = info.py signon.py
PROJ_CPPTESTS = info inject many
PROJ_CTESTS = c_info c_many
PROJ_TESTS = $(PROJ_PYTESTS) $(PROJ_CPPTESTS) $(PROJ_CTESTS)
PROJ_LIBRARIES = -ljournal

ifeq (AIX, ${findstring AIX, $(PLATFORM_ID)}) 
    PROJ_LCXX_FLAGS += -brtl
endif

#--------------------------------------------------------------------------
#

all: $(PROJ_TESTS)

test:
	for test in $(PROJ_TESTS) ; do ./$${test}; done

release: tidy
	cvs release .

update: clean
	cvs update .

#--------------------------------------------------------------------------
#

info: info.cc $(BLD_LIBDIR)/libjournal.$(EXT_SO)
	$(CXX) $(CXXFLAGS) $(LCXXFLAGS) -o $@ info.cc $(PROJ_LIBRARIES)

inject: inject.cc $(BLD_LIBDIR)/libjournal.$(EXT_SO)
	$(CXX) $(CXXFLAGS) $(LCXXFLAGS) -o $@ inject.cc $(PROJ_LIBRARIES)

many: many.cc $(BLD_LIBDIR)/libjournal.$(EXT_SO)
	$(CXX) $(CXXFLAGS) $(LCXXFLAGS) -o $@ many.cc $(PROJ_LIBRARIES)

c_info: c_info.c $(BLD_LIBDIR)/libjournal.$(EXT_SO)
	$(CC) $(CFLAGS) -c c_info.c
	$(CXX) $(LCXXFLAGS) -o $@ c_info.o $(PROJ_LIBRARIES)
	$(RM) $(RMFLAGS) c_info.o

c_many: c_many.c $(BLD_LIBDIR)/libjournal.$(EXT_SO)
	$(CC) $(CFLAGS) -c c_many.c
	$(CXX) $(LCXXFLAGS) -o $@ c_many.o $(PROJ_LIBRARIES)
	$(RM) $(RMFLAGS) c_many.o

# version
# $Id: Make.mm,v 1.7 2003/05/07 20:30:20 cummings Exp $

# End of file
