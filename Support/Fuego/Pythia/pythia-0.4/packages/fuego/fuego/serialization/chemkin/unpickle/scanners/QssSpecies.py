#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from BaseScanner import BaseScanner


class QssSpecies(BaseScanner):


    def _tokenClasses(self):
        from fuego.serialization.chemkin.unpickle.tokens import qss_speciesTokenClasses
        return qss_speciesTokenClasses()


# version
__id__ = "$Id$"

#
# End of file
