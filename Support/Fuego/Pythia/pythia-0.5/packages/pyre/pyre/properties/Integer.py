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

from .Property import Property


class Integer(Property):
    def __init__(
        self, name, default=0, public=None, validator=None, tip="", doc=""
    ):
        super(Integer, self).__init__(
            name, default, public, validator, tip, doc
        )
        self.type = "integer"
        return

    def _cast(self, value):
        return int(value)


# version
__id__ = "$Id$"

# End of file
