#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from __future__ import absolute_import

from .BaseScanner import BaseScanner


class Reactions(BaseScanner):
    def _tokenClasses(self):
        from fuego.serialization.chemkin.unpickle.tokens import (
            reactionTokenClasses,
        )

        return reactionTokenClasses()


# version
__id__ = "$Id$"

#
# End of file
