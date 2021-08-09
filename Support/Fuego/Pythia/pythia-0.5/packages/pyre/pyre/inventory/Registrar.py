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
from pyre.properties.Property import Property


class Registrar(type):
    def __init__(cls, name, bases, dict):
        super(Registrar, cls).__init__(name, bases, dict)

        propertyRegistry = {}
        facilityRegistry = {}
        componentRegistry = {}

        # add to the registry properties inherited from ancestors
        bases = list(cls.__bases__)
        bases.reverse()
        for base in bases:
            try:
                propertyRegistry.update(base._propertyRegistry)
            except AttributeError:
                pass

            try:
                facilityRegistry.update(base._facilityRegistry)
            except AttributeError:
                pass

            try:
                componentRegistry.update(base._componentRegistry)
            except AttributeError:
                pass

        # register the traits defined by this class
        for trait in cls.inventory:
            name = trait.name
            public = trait.public

            if isinstance(trait, Property):
                propertyRegistry[public] = trait
            elif isinstance(trait, Facility):
                facilityRegistry[public] = trait
                # if trait.default is not None:
                # componentRegistry[trait.default.name] = trait
                # else:
                # componentRegistry[public] = trait

            setattr(cls, name, trait)

        # update the class record
        cls.inventory = ()
        cls._propertyRegistry = propertyRegistry
        cls._facilityRegistry = facilityRegistry
        cls._componentRegistry = componentRegistry

        return


# version
__id__ = "$Id$"

# End of file
