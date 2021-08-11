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

from __future__ import absolute_import

from pyre.xml.DocumentNode import DocumentNode


class Configuration(DocumentNode):
    def options(self):
        return self._options

    def onWeaver(self, options):
        self._options = options
        return

    def __init__(self, source):
        DocumentNode.__init__(self, source, documentNodes())
        return


# helpers
def documentNodes():

    from .Author import Author
    from .BannerCharacter import BannerCharacter
    from .BannerWidth import BannerWidth
    from .Copyright import Copyright
    from .Mill import Mill
    from .Organization import Organization
    from .Weaver import Weaver

    nodes = [
        Weaver,
        Mill,
        Author,
        Organization,
        Copyright,
        BannerCharacter,
        BannerWidth,
    ]
    return nodes


# version
__id__ = "$Id$"

#  End of file
