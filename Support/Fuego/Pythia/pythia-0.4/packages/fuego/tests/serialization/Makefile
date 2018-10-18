# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2001  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT =

#--------------------------------------------------------------------------
#

all: hoxy


.PHONY: cmech
cmech: grimech.c
	$(CC) -O2 grimech.c -o grimech -lm
	grimech


.PHONY: grimech.c
grimech.c:
	pickle.py --file=grimech.c


.PHONY: hoxy
hoxy: hoxy.c
	$(CC) -O2 hoxy.c -o hoxy -lm
	hoxy


.PHONY: hoxy.c
hoxy.c:
	pickle.py --file=hoxy.c --mechanism=HydrogenOxygen.ck2

#--------------------------------------------------------------------------
#

clean::
	@-rm hoxy* grimech* core
	@find * -name \*.bak -exec rm {} \;
	@find * -name \*.pyc -exec rm {} \;
	@find * -name \*~ -exec rm {} \;

# version
# $Id$

#
# End of file
