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


from weaver.mills.XMLMill import XMLMill


class Renderer(XMLMill):


    def onRegistry(self, registry):

        # bail out of empty registries
        if not registry.properties and not registry.facilities:
            return
        
        self._indent()
        self._write('<facility name="%s">' % registry.name)

        self._indent()
        for property, value in registry.properties.iteritems():
            self._write('<property name="%s" value="%s"/>' % (property, value))
        self._outdent()

        for facility in registry.facilities.itervalues():
            self.onRegistry(facility)

        self._write('</facility>')
        self._outdent()

        return


    def render(self, registry):
        document = self.pickle(registry)
        return document


    def __init__(self):
        XMLMill.__init__(self)
        return


    def _renderDocument(self, registry):
        self._rep += ['', '<!DOCTYPE registry>', '', '<registry>', '' ]
        self.onRegistry(registry)
        self._rep += ['', '</registry>']

        return
    

# version
__id__ = "$Id$"

# End of file 
