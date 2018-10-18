#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


from Device import Device


class File(Device):


    def __init__(self, logfile):
        Device.__init__(self)
        self.file = logfile
        return


    def _write(self, message):
        for line in message:
            self.file.write("%s\n" % line)

        self.file.flush()

        return


# version
__id__ = "$Id$"

#  End of file 
