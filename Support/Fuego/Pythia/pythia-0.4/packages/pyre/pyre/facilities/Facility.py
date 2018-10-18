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


from pyre.inventory.Trait import Trait


class Facility(Trait):


    def __init__(
        self, name, default=None, public=None, 
        binder=None, requirements=None, tip="", doc="" ):

        Trait.__init__(self, name, default=default, public=public, tip=tip, doc=doc)
        self.type = "facility"

        if binder is None:
            from ScriptBinder import ScriptBinder
            binder = ScriptBinder()
        self.binder = binder
        
        if requirements is None:
            requirements = []
        self.requirements = requirements
            
        return


    def __set__(self, instance, component):

        if isinstance(component, str):
            component = self.binder.bind(self.name, component)

        # binding succeeded; install the component as the value of the facility
        instance.__dict__[self.mangled] = component

        return


# version
__id__ = "$Id$"

# End of file 
