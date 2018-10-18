#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from pyre.facilities.Facility import Facility


class JournalFacility(Facility):


    def __init__(self, default=None):
        if default is None:
            from Journal import Journal
            default = Journal()
            
        Facility.__init__(self, name="journal", default=default)
        return


# version
__id__ = "$Id$"

# End of file 
