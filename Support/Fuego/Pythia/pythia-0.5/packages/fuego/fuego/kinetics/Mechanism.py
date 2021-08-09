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

from builtins import object


class Mechanism(object):
    def reaction(self, id=None):
        return self._reactions.find(id)

    def __init__(self, mixture, mechanism):
        from .ReactionSet import ReactionSet

        self._reactions = ReactionSet(mixture, mechanism.reaction())
        return


# version
__id__ = "$Id$"

# End of file
