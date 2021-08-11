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

from .AbstractNode import AbstractNode


class Organization(AbstractNode):

    tag = "organization"

    def content(self, text):
        self._organization = text
        return

    def notify(self, parent):
        parent.onOrganization(self._organization)
        return

    def __init__(self, root, attributes):
        AbstractNode.__init__(self, root, attributes)
        self._organization = ""
        return


# version
__id__ = "$Id$"

#  End of file
