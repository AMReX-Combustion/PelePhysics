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

include chemkin/default.def

TESTS = productionRates progressRates

all: $(TESTS)

CLEAN += $(TESTS)

LIBS = $(EXTERNAL_LIBS)


productionRates: productionRates.f
	$(F77) $(F77FLAGS) $(LF77FLAGS) productionRates.f -o productionRates $(LIBS)


progressRates: progressRates.f
	$(F77) $(F77FLAGS) $(LF77FLAGS) progressRates.f -o progressRates $(LIBS)

#
# End of file
