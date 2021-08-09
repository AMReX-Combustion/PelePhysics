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

from pyre.applications.Registrar import Registrar as BaseRegistrar


class Registrar(BaseRegistrar):


    def retrieve(self, name):
        return self.registry.get(name)


    def register(self, entity, id, aliases=[]):
        self.registered.append(id)
        self.registry[id] = entity
        for alias in aliases:
            self.registry[alias] = entity
        return


    def init(self):
        BaseRegistrar.init(self)
        self.registered = []
        return


# version
__id__ = "$Id$"

#  End of file 
