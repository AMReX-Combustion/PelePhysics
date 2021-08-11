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

import pyre.util.bool

from .Property import Property


class Bool(Property):
    def __init__(
        self, name, default=False, public=None, validator=None, tip="", doc=""
    ):
        super(Bool, self).__init__(name, default, public, validator, tip, doc)
        self.type = "bool"
        return

    def _cast(self, value):
        if isinstance(value, str):
            return pyre.util.bool.bool(value)

        return bool(value)


# version
__id__ = "$Id$"

# End of file
