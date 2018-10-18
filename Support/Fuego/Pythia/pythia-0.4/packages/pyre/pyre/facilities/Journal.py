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


from Facility import Facility
from ScriptBinder import ScriptBinder

class Journal(Facility, ScriptBinder):


    def __init__(self, default=None):
        if default is None:
            import journal
            default = journal.journal()
            
        Facility.__init__(self, "journal", default=default)
        ScriptBinder.__init__(self)
        return


# version
__id__ = "$Id$"

# End of file 
