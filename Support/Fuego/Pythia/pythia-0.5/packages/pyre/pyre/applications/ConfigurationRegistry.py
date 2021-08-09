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


from pyre.inventory.Registry import Registry as BaseRegistry


class ConfigurationRegistry(BaseRegistry):


    def __init__(self, name):
        BaseRegistry.__init__(self, name)

        self.help = False
        self.orphans = []
        self.unknown = {}

        return


# version
__id__ = "$Id$"

#  End of file 
