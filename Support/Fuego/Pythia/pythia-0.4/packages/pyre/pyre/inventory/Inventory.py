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


class Inventory(object):


    # class data
    delegate = None
    inventory = ()


    def defaults(self):
        for facility in self._facilityRegistry:
            component = self.__getattribute__(facility)
            component.defaults()
        return


    def init(self, parent):
        for facility in self._facilityRegistry:
            component = self.__getattribute__(facility)
            component.init(parent)
        return


    def configure(self, registry):
        unknown = self.configureProperties(registry)
        unknown += self.configureComponents(registry)

        report = [ ("%s.%s" % (registry.name, name), value) for name, value in unknown ]

        return report


    def fini(self):
        for facility in self._facilityRegistry:
            component = self.__getattribute__(facility)
            component.fini()
        return


    def configureProperties(self, registry):
        unknown = []

        for name, value in registry.properties.iteritems():
            # attempt to look up the public name in the property registries
            try:
                prop = self._propertyRegistry[name]
            except KeyError:
                try:
                    prop = self._facilityRegistry[name]
                except KeyError:
                    unknown.append((name, value))
                    continue

            prop.__set__(self, value)

        return unknown         
    

    def configureComponents(self, registry):
        unknown = []

        # build a dictionary of all known facilities and components
        table = {}
        for name in self._facilityRegistry:
            component = self.__getattribute__(name)
            table[name] = component
            table[component.name] = component

        # iterate over the registry entries
        for name, conf in registry.facilities.iteritems():
            try:
                component = table[name]
            except KeyError:
                unknown += conf.list()
                continue

            unknown += component.configure(conf)
        
        return unknown


    def facilities(self):
        return dict(self._facilityRegistry)


    def properties(self):
        return dict(self._propertyRegistry)


    def components(self):
        table = {}
        for name, facility in self._facilityRegistry.iteritems():
            table[name] = self.__getattribute__(facility.name)

        return table


    def state(self, registry):
        for property in self._propertyRegistry.itervalues():
            registry.properties[property.public] = self.__getattribute__(property.name)

        for facility in self._facilityRegistry.itervalues():
            component = self.__getattribute__(facility.name)
            node = registry.open(component.name)
            component.state(node)
            
        return registry


    # delegation implementation
    def __getattribute__(self, attribute):
        # attempt to get the requested attribute
        try:
            return super(Inventory, self).__getattribute__(attribute)
        except AttributeError:
            pass

        # if we don't have one, see whether there is someone to whom we can delegate this request
        try:
            return self.delegate.__getattribute__(attribute)
        except AttributeError:
            pass

        # raise an exception
        raise AttributeError, "'%s' object has no property named '%s'" % (
            self.__class__.__name__, attribute)

    # implementation details: data
    from Registrar import Registrar
    __metaclass__ = Registrar


# version
__id__ = "$Id$"

# End of file 
