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


from __future__ import absolute_import
from pyre.xml.DocumentNode import DocumentNode


class Document(DocumentNode):


    def onRegistry(self, registry):
        self.registry = registry
        return


    def __init__(self, source):
        DocumentNode.__init__(self, source, documentNodes())
        self.registry = None
        return


# helpers

def documentNodes():

    from .Registry import Registry
    from .Facility import Facility
    from .Property import Property

    nodes = [ Registry, Facility, Property ]

    return nodes

    
# version
__id__ = "$Id$"

# End of file 
