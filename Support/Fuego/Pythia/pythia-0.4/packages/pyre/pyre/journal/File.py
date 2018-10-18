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


from Device import Device


class File(Device):


    def createDevice(self):
        from journal.File import File
        logfile = file(self.inventory.name, "w")
        return File(logfile)


    def __init__(self):
        Device.__init__(self, "file")
        return


    class Inventory(Device.Inventory):


        import pyre.properties


        inventory = (
            pyre.properties.str("name", default="journal.log"),
            )


# version
__id__ = "$Id$"

# End of file 
